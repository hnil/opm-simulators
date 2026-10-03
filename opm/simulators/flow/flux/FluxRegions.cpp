/*
  Copyright 2026 Equinor ASA

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/simulators/flow/flux/FluxRegions.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <set>
#include <unordered_set>

namespace Opm {

namespace {

struct GlobalIJK {
    int i = 0;
    int j = 0;
    int k = 0;
};

GlobalIJK decodeGlobal(const std::array<int, 3>& dims, const int globalIndex)
{
    const int k = globalIndex / (dims[0] * dims[1]);
    const int rem = globalIndex % (dims[0] * dims[1]);
    const int j = rem / dims[0];
    const int i = rem % dims[0];
    return {i + 1, j + 1, k + 1};
}

int encodeGlobal(const std::array<int, 3>& dims, const GlobalIJK& ijk)
{
    return FluxRegions::cartesianIndex(dims, ijk.i - 1, ijk.j - 1, ijk.k - 1);
}

bool inBounds(const std::array<int, 3>& dims, const GlobalIJK& ijk)
{
    return (ijk.i >= 1 && ijk.i <= dims[0])
        && (ijk.j >= 1 && ijk.j <= dims[1])
        && (ijk.k >= 1 && ijk.k <= dims[2]);
}

GlobalIJK shifted(const GlobalIJK& ijk, FaceDir::DirEnum dir)
{
    switch (dir) {
    case FaceDir::XPlus:
        return {ijk.i + 1, ijk.j, ijk.k};
    case FaceDir::XMinus:
        return {ijk.i - 1, ijk.j, ijk.k};
    case FaceDir::YPlus:
        return {ijk.i, ijk.j + 1, ijk.k};
    case FaceDir::YMinus:
        return {ijk.i, ijk.j - 1, ijk.k};
    case FaceDir::ZPlus:
        return {ijk.i, ijk.j, ijk.k + 1};
    case FaceDir::ZMinus:
        return {ijk.i, ijk.j, ijk.k - 1};
    case FaceDir::Unknown:
        break;
    }

    return ijk;
}

FaceDir::DirEnum directionToNeighbour(const GlobalIJK& from, const GlobalIJK& to)
{
    if (to.i == from.i + 1 && to.j == from.j && to.k == from.k) {
        return FaceDir::XPlus;
    }
    if (to.i == from.i - 1 && to.j == from.j && to.k == from.k) {
        return FaceDir::XMinus;
    }
    if (to.j == from.j + 1 && to.i == from.i && to.k == from.k) {
        return FaceDir::YPlus;
    }
    if (to.j == from.j - 1 && to.i == from.i && to.k == from.k) {
        return FaceDir::YMinus;
    }
    if (to.k == from.k + 1 && to.i == from.i && to.j == from.j) {
        return FaceDir::ZPlus;
    }
    if (to.k == from.k - 1 && to.i == from.i && to.j == from.j) {
        return FaceDir::ZMinus;
    }

    return FaceDir::Unknown;
}

std::vector<FaceDir::DirEnum> cartesianDirections()
{
    return {
        FaceDir::XMinus,
        FaceDir::XPlus,
        FaceDir::YMinus,
        FaceDir::YPlus,
        FaceDir::ZMinus,
        FaceDir::ZPlus,
    };
}

std::unordered_set<int> makeSelectionSet(const FluxRegions::Region& region)
{
    return {region.selectedGlobalCells.begin(), region.selectedGlobalCells.end()};
}

std::map<int, std::vector<int>> buildGlobalToLocal(const FluxRegions::Region& region)
{
    std::map<int, std::vector<int>> map;
    for (std::size_t local = 0; local < region.localToGlobal.size(); ++local) {
        const auto global = region.localToGlobal[local];
        if (global >= 0) {
            map[global].push_back(static_cast<int>(local));
        }
    }
    return map;
}

void fillCartesianBoundaryFaces(const std::array<int, 3>& dims,
                                FluxRegions::Region& region,
                                const std::unordered_set<int>& selected,
                                const std::vector<int>& actnum)
{
    const auto isActive = [&actnum](const int global)
    {
        return actnum.empty() || (actnum[global] != 0);
    };

    const auto globalToLocal = buildGlobalToLocal(region);
    for (const auto interiorGlobal : region.selectedGlobalCells) {
        const auto interiorIt = globalToLocal.find(interiorGlobal);
        if (interiorIt == globalToLocal.end() || interiorIt->second.empty()) {
            continue;
        }

        const int interiorLocal = interiorIt->second.front();
        const auto interiorIJK = decodeGlobal(dims, interiorGlobal);

        for (const auto dir : cartesianDirections()) {
            const auto nIJK = shifted(interiorIJK, dir);
            if (!inBounds(dims, nIJK)) {
                continue;
            }

            const auto neighbourGlobal = encodeGlobal(dims, nIJK);
            if (selected.count(neighbourGlobal) != 0) {
                continue;
            }

            // Nothing flows across a face towards an inactive cell, so it is a
            // no-flow boundary rather than one the parent run has to supply a
            // rate or a pressure for.
            if (!isActive(neighbourGlobal)) {
                continue;
            }

            region.boundaryFaces.push_back(
                FluxRegions::BoundaryFace{
                    interiorLocal,
                    interiorGlobal,
                    neighbourGlobal,
                    dir,
                    false,
                });
        }
    }
}

void fillNncBoundaryFaces(const std::array<int, 3>& dims,
                          FluxRegions::Region& region,
                          const std::unordered_set<int>& selected,
                          const std::vector<std::array<int, 2>>& nncConnections)
{
    const auto globalToLocal = buildGlobalToLocal(region);
    std::set<std::tuple<int, int, int>> seen;

    auto addFace = [&](const int interiorGlobal, const int exteriorGlobal) {
        if (selected.count(interiorGlobal) == 0 || selected.count(exteriorGlobal) != 0) {
            return;
        }

        const auto localIt = globalToLocal.find(interiorGlobal);
        if (localIt == globalToLocal.end() || localIt->second.empty()) {
            return;
        }

        const auto interiorIJK = decodeGlobal(dims, interiorGlobal);
        const auto exteriorIJK = decodeGlobal(dims, exteriorGlobal);
        auto dir = directionToNeighbour(interiorIJK, exteriorIJK);

        const auto key = std::make_tuple(localIt->second.front(), exteriorGlobal, static_cast<int>(dir));
        if (!seen.insert(key).second) {
            return;
        }

        region.boundaryFaces.push_back(
            FluxRegions::BoundaryFace{
                localIt->second.front(),
                interiorGlobal,
                exteriorGlobal,
                dir,
                true,
            });
    };

    for (const auto& pair : nncConnections) {
        addFace(pair[0], pair[1]);
        addFace(pair[1], pair[0]);
    }
}

void sortBoundaryFaces(std::vector<FluxRegions::BoundaryFace>& faces)
{
    std::sort(faces.begin(), faces.end(),
              [](const auto& lhs, const auto& rhs) {
                  if (lhs.interiorGlobalCell != rhs.interiorGlobalCell) {
                      return lhs.interiorGlobalCell < rhs.interiorGlobalCell;
                  }
                  if (lhs.isNnc != rhs.isNnc) {
                      return lhs.isNnc < rhs.isNnc;
                  }
                  if (lhs.direction != rhs.direction) {
                      return lhs.direction < rhs.direction;
                  }
                  return lhs.exteriorGlobalCell < rhs.exteriorGlobalCell;
              });
}

void buildBoundaryFaces(const std::array<int, 3>& dims,
                        FluxRegions::Region& region,
                        const std::vector<int>& actnum,
                        const std::vector<std::array<int, 2>>& nncConnections)
{
    const auto selected = makeSelectionSet(region);
    fillCartesianBoundaryFaces(dims, region, selected, actnum);
    fillNncBoundaryFaces(dims, region, selected, nncConnections);
    sortBoundaryFaces(region.boundaryFaces);
}

} // namespace

int FluxRegions::cartesianIndex(const std::array<int, 3>& dims,
                                const int i,
                                const int j,
                                const int k)
{
    return i + dims[0] * (j + dims[1] * k);
}

int FluxRegions::localBoxIndex(const Box& box, const int i, const int j, const int k)
{
    return (i - box.i1) + box.nx() * ((j - box.j1) + box.ny() * (k - box.k1));
}

std::vector<FluxRegions::Region> FluxRegions::extract(const std::array<int, 3>& dims,
                                                      const std::vector<int>& regionValues)
{
    return extract(dims, regionValues, {}, {});
}

std::vector<FluxRegions::Region> FluxRegions::extract(const std::array<int, 3>& dims,
                                                      const std::vector<int>& regionValues,
                                                      const std::vector<std::array<int, 2>>& nncConnections)
{
    return extract(dims, regionValues, {}, nncConnections);
}

std::vector<FluxRegions::Region> FluxRegions::extract(const std::array<int, 3>& dims,
                                                      const std::vector<int>& regionValues,
                                                      const std::vector<int>& actnum,
                                                      const std::vector<std::array<int, 2>>& nncConnections)
{
    const auto numCells = dims[0] * dims[1] * dims[2];
    if (static_cast<int>(regionValues.size()) != numCells) {
        OPM_THROW(std::invalid_argument,
                  "FluxRegions::extract(): regionValues size must match grid dimensions");
    }

    if (!actnum.empty() && (static_cast<int>(actnum.size()) != numCells)) {
        OPM_THROW(std::invalid_argument,
                  "FluxRegions::extract(): actnum size must match grid dimensions");
    }

    // The region map is dimensioned over every cell of the parent grid, so it
    // also assigns a region to cells that are inactive. Those cells hold no
    // fluid and take no part in the flow, so they must not become part of the
    // region: the sector's active cells are the region combined with ACTNUM,
    // and the local-to-active mapping the consumer builds from this file
    // assumes exactly that set.
    const auto isActive = [&actnum](const int global)
    {
        return actnum.empty() || (actnum[global] != 0);
    };

    struct PendingRegion {
        Box box;
        std::vector<int> globalCells;
    };

    std::map<int, PendingRegion> pending;

    for (int k = 0; k < dims[2]; ++k) {
        for (int j = 0; j < dims[1]; ++j) {
            for (int i = 0; i < dims[0]; ++i) {
                const auto globalIndex = cartesianIndex(dims, i, j, k);
                const auto regionId = regionValues[globalIndex];
                if (regionId <= 0 || !isActive(globalIndex)) {
                    continue;
                }

                auto& region = pending[regionId];
                if (region.globalCells.empty()) {
                    region.box = Box{dims[0] + 1, 0, dims[1] + 1, 0, dims[2] + 1, 0};
                }
                region.box.i1 = std::min(region.box.i1, i + 1);
                region.box.i2 = std::max(region.box.i2, i + 1);
                region.box.j1 = std::min(region.box.j1, j + 1);
                region.box.j2 = std::max(region.box.j2, j + 1);
                region.box.k1 = std::min(region.box.k1, k + 1);
                region.box.k2 = std::max(region.box.k2, k + 1);
                region.globalCells.push_back(globalIndex);
            }
        }
    }

    std::vector<Region> regions;
    regions.reserve(pending.size());

    for (auto& [regionId, pendingRegion] : pending) {
        Region region;
        region.regionId = regionId;
        region.box = pendingRegion.box;
        region.selectedGlobalCells = std::move(pendingRegion.globalCells);

        region.localToGlobal.assign(region.box.nx() * region.box.ny() * region.box.nz(), -1);
        for (const auto globalIndex : region.selectedGlobalCells) {
            const int k = globalIndex / (dims[0] * dims[1]);
            const int rem = globalIndex % (dims[0] * dims[1]);
            const int j = rem / dims[0];
            const int i = rem % dims[0];
            region.localToGlobal[localBoxIndex(region.box, i + 1, j + 1, k + 1)] = globalIndex;
        }

        buildBoundaryFaces(dims, region, actnum, nncConnections);

        regions.push_back(std::move(region));
    }

    return regions;
}

int FluxRegions::uniqueSelectedRegion(const std::vector<int>& regionValues)
{
    std::set<int> regionIds;
    for (const auto value : regionValues) {
        if (value > 0) {
            regionIds.insert(value);
        }
    }

    if (regionIds.empty()) {
        OPM_THROW(std::invalid_argument,
                  "FluxRegions::uniqueSelectedRegion(): FLUXNUM must contain at least one positive region id");
    }

    if (regionIds.size() != 1U) {
        OPM_THROW(std::invalid_argument,
                  "FluxRegions::uniqueSelectedRegion(): FLUXNUM must contain exactly one positive region id");
    }

    return *regionIds.begin();
}

std::vector<int> FluxRegions::buildActnum(const std::vector<int>& regionValues,
                                          const int regionId)
{
    std::vector<int> actnum(regionValues.size(), 0);
    for (std::size_t index = 0; index < regionValues.size(); ++index) {
        actnum[index] = (regionValues[index] == regionId) ? 1 : 0;
    }

    return actnum;
}

std::vector<int> FluxRegions::buildActnum(const std::vector<int>& regionValues,
                                          const std::vector<int>& regionIds)
{
    std::set<int> selectedIds;
    for (const auto regionId : regionIds) {
        if (regionId > 0) {
            selectedIds.insert(regionId);
        }
    }

    if (selectedIds.empty()) {
        OPM_THROW(std::invalid_argument,
                  "FluxRegions::buildActnum(): FLUXREG must list at least one positive region id");
    }

    std::set<int> availableIds;
    for (const auto value : regionValues) {
        if (value > 0) {
            availableIds.insert(value);
        }
    }

    for (const auto regionId : selectedIds) {
        if (!availableIds.count(regionId)) {
            OPM_THROW(std::invalid_argument,
                      "FluxRegions::buildActnum(): FLUXREG selects a region id not present in FLUXNUM");
        }
    }

    std::vector<int> actnum(regionValues.size(), 0);
    for (std::size_t index = 0; index < regionValues.size(); ++index) {
        actnum[index] = selectedIds.count(regionValues[index]) ? 1 : 0;
    }

    return actnum;
}

} // namespace Opm