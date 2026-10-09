/*
  Copyright 2026 Equinor ASA.

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
#ifndef OPM_WELL_ZONE_LGR_HPP
#define OPM_WELL_ZONE_LGR_HPP

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm::Parameters {

//! \brief Nested refinement rings around wells (--well-refine).
struct WellRefine { static constexpr auto value = ""; };

} // namespace Opm::Parameters

namespace Opm {

//! One well-zone spec: which wells, the per-ring factor, the ring widths
//! (outermost first, in coarse cells beyond the next ring), and the K extent.
struct WellZoneSpec
{
    std::string pattern;
    std::array<int,3> factor{};
    std::vector<int> rings;         //!< outermost first
    bool allLayers = false;         //!< K = whole grid instead of the perforated layers
    int extraLayers = 0;            //!< layers added above and below the perforated ones
};

/// Parse "--well-refine". Zones separated by ';'; each is
///   PATTERN:FX,FY,FZ:rings=W1[,W2,...][:layers=all][:klayers=N]
/// e.g. "INJ*:3,3,1:rings=2,1" refines a 3-cell ring by 3x3 in I,J and, nested
/// inside it, a 1-cell ring by another 3x3, so the well column ends up 9x9.
inline std::vector<WellZoneSpec> parseWellZoneSpec(const std::string& spec)
{
    std::vector<WellZoneSpec> zones;
    std::stringstream zoneStream(spec);
    std::string one;
    while (std::getline(zoneStream, one, ';')) {
        if (one.find_first_not_of(" \t") == std::string::npos) {
            continue;
        }
        std::vector<std::string> parts;
        std::stringstream ps(one);
        std::string p;
        while (std::getline(ps, p, ':')) {
            parts.push_back(p);
        }
        auto bad = [&one](const std::string& why) {
            throw std::invalid_argument(
                "--well-refine: zone '" + one + "': " + why +
                ". Expected PATTERN:FX,FY,FZ:rings=W1[,W2,...][:layers=all][:klayers=N].");
        };
        if (parts.size() < 3) {
            bad("too few fields");
        }
        WellZoneSpec z;
        z.pattern = parts[0];
        {
            std::stringstream fs(parts[1]);
            char c1 = 0, c2 = 0;
            if (!(fs >> z.factor[0] >> c1 >> z.factor[1] >> c2 >> z.factor[2]) || c1 != ',' || c2 != ',') {
                bad("factor must be FX,FY,FZ");
            }
            for (int f : z.factor) {
                if (f <= 0) {
                    bad("factors must be positive");
                }
            }
        }
        for (std::size_t i = 2; i < parts.size(); ++i) {
            const auto eq = parts[i].find('=');
            const auto key = parts[i].substr(0, eq);
            const auto val = (eq == std::string::npos) ? std::string{} : parts[i].substr(eq + 1);
            if (key == "rings") {
                std::stringstream rs(val);
                std::string w;
                while (std::getline(rs, w, ',')) {
                    const int width = std::stoi(w);
                    if (width < 1) {
                        bad("ring widths must be at least 1 coarse cell");
                    }
                    z.rings.push_back(width);
                }
            }
            else if (key == "layers") {
                if (val != "all" && val != "perforated") {
                    bad("layers must be 'all' or 'perforated'");
                }
                z.allLayers = (val == "all");
            }
            else if (key == "klayers") {
                z.extraLayers = std::stoi(val);
            }
            else {
                bad("unknown field '" + key + "'");
            }
        }
        if (z.rings.empty()) {
            bad("at least one ring is required");
        }
        zones.push_back(std::move(z));
    }
    return zones;
}

namespace detail {

struct CellBox
{
    std::array<int,3> lo{ std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max() };
    std::array<int,3> hi{ std::numeric_limits<int>::min(), std::numeric_limits<int>::min(), std::numeric_limits<int>::min() };
    bool empty() const { return lo[0] > hi[0]; }
    void add(const std::array<int,3>& c)
    {
        for (int d = 0; d < 3; ++d) {
            lo[d] = std::min(lo[d], c[d]);
            hi[d] = std::max(hi[d], c[d]);
        }
    }
    bool overlapsOrTouches(const CellBox& o) const
    {
        for (int d = 0; d < 3; ++d) {
            if (hi[d] + 1 < o.lo[d] || o.hi[d] + 1 < lo[d]) {
                return false;
            }
        }
        return true;
    }
    void merge(const CellBox& o)
    {
        for (int d = 0; d < 3; ++d) {
            lo[d] = std::min(lo[d], o.lo[d]);
            hi[d] = std::max(hi[d], o.hi[d]);
        }
    }
};

//! Merge boxes that overlap or touch, repeatedly, so no two survivors do.
inline void mergeBoxes(std::vector<CellBox>& boxes)
{
    bool merged = true;
    while (merged) {
        merged = false;
        for (std::size_t a = 0; a < boxes.size() && !merged; ++a) {
            for (std::size_t b = a + 1; b < boxes.size(); ++b) {
                if (boxes[a].overlapsOrTouches(boxes[b])) {
                    boxes[a].merge(boxes[b]);
                    boxes.erase(boxes.begin() + b);
                    merged = true;
                    break;
                }
            }
        }
    }
}

} // namespace detail

/// Build the nested ring requests for every zone.  The perforated cells of all
/// wells matching a zone's pattern, over the whole schedule, are the seed; the
/// outermost ring is the seed dilated by the sum of the widths laterally and by
/// one layer per inner ring in K (the innermost covers the perforated layers,
/// or every ring covers every layer with layers=all), each inner ring by one
/// width and one layer less, nested inside the previous one with the same factor.  Boxes
/// of one zone that overlap or touch are merged, and a ring that would touch
/// its parent's boundary after clamping to the grid is dropped with a warning.
inline std::vector<Opm::Refinement::BlockRefinement>
wellZoneRefinements(const std::vector<WellZoneSpec>& zones,
                    const Schedule& schedule,
                    const std::array<int,3>& cartDims)
{
    using Opm::Refinement::BlockRefinement;
    std::vector<BlockRefinement> requests;
    int zoneIdx = 0;

    for (const auto& zone : zones) {
        ++zoneIdx;
        // Seed boxes: one per well, from every report step's connections.
        std::vector<detail::CellBox> seeds;
        std::vector<std::string> seen;
        for (std::size_t step = 0; step < schedule.size(); ++step) {
            for (const auto& name : schedule.wellNames(zone.pattern, step)) {
                const auto& well = schedule.getWell(name, step);
                auto it = std::find(seen.begin(), seen.end(), name);
                if (it == seen.end()) {
                    seen.push_back(name);
                    seeds.emplace_back();
                    it = std::prev(seen.end());
                }
                auto& box = seeds[std::distance(seen.begin(), it)];
                for (const auto& conn : well.getConnections()) {
                    box.add({ conn.getI(), conn.getJ(), conn.getK() });
                }
            }
        }
        seeds.erase(std::remove_if(seeds.begin(), seeds.end(),
                                   [](const detail::CellBox& b) { return b.empty(); }),
                    seeds.end());
        if (seeds.empty()) {
            OpmLog::warning(fmt::format("--well-refine: no perforated well matches '{}'; "
                                        "zone skipped.", zone.pattern));
            continue;
        }

        const int nRings = static_cast<int>(zone.rings.size());
        // Ring r (0 = outermost) extends beyond the seed by the widths of rings r..end.
        std::vector<int> reach(nRings, 0);
        for (int r = nRings - 1; r >= 0; --r) {
            reach[r] = zone.rings[r] + ((r + 1 < nRings) ? reach[r + 1] : 0);
        }

        // Coarse-grid boxes per ring, merged per ring.  Every ring dilates the
        // same seeds, so a merged outer ring contains the merged inner rings.
        std::vector<std::vector<detail::CellBox>> rings(nRings);
        for (int r = 0; r < nRings; ++r) {
            for (const auto& seed : seeds) {
                detail::CellBox b = seed;
                for (int d = 0; d < 2; ++d) {
                    b.lo[d] = std::max(0, seed.lo[d] - reach[r]);
                    b.hi[d] = std::min(cartDims[d] - 1, seed.hi[d] + reach[r]);
                }
                if (zone.allLayers) {
                    b.lo[2] = 0;
                    b.hi[2] = cartDims[2] - 1;
                }
                else {
                    // A nested ring must be interior to its parent in K too, so
                    // each ring reaches one layer beyond the ring inside it; the
                    // innermost covers exactly the perforated layers (plus any
                    // extra asked for).
                    const int kReach = zone.extraLayers + (nRings - 1 - r);
                    b.lo[2] = std::max(0, seed.lo[2] - kReach);
                    b.hi[2] = std::min(cartDims[2] - 1, seed.hi[2] + kReach);
                }
                rings[r].push_back(b);
            }
            detail::mergeBoxes(rings[r]);
        }

        // Emit: outer ring boxes on GLOBAL, each inner box nested in the outer
        // box containing it, in that box's refined coordinates (a coarse cell
        // is unitsPerCoarse refined cells of the parent).
        struct Emitted { detail::CellBox coarse; std::string name; std::array<int,3> unitsPerCoarse; };
        std::vector<Emitted> previous;
        for (int r = 0; r < nRings; ++r) {
            std::vector<Emitted> current;
            int boxIdx = 0;
            for (const auto& box : rings[r]) {
                ++boxIdx;
                BlockRefinement req;
                req.name = fmt::format("WZ{}R{}B{}", zoneIdx, r + 1, boxIdx);
                req.cellsPerDim = zone.factor;
                if (r == 0) {
                    for (int d = 0; d < 3; ++d) {
                        req.startIJK[d] = box.lo[d];
                        req.endIJK[d] = box.hi[d] + 1;
                    }
                    current.push_back({ box, req.name, zone.factor });
                    requests.push_back(std::move(req));
                    continue;
                }
                const Emitted* parent = nullptr;
                for (const auto& p : previous) {
                    bool inside = true;
                    for (int d = 0; d < 3; ++d) {
                        inside = inside && box.lo[d] >= p.coarse.lo[d] && box.hi[d] <= p.coarse.hi[d];
                    }
                    if (inside) {
                        parent = &p;
                        break;
                    }
                }
                if (parent == nullptr) {
                    throw std::logic_error("--well-refine: an inner ring is not contained "
                                           "in any outer ring box");
                }
                // Strictly interior to the parent: the builder needs at least one
                // refined parent cell on every side, which clamping at the grid
                // boundary can take away; shrink the ring there.
                detail::CellBox inner = box;
                bool shrunk = false;
                for (int d = 0; d < 3; ++d) {
                    if (inner.lo[d] <= parent->coarse.lo[d]) {
                        inner.lo[d] = parent->coarse.lo[d] + 1;
                        shrunk = true;
                    }
                    if (inner.hi[d] >= parent->coarse.hi[d]) {
                        inner.hi[d] = parent->coarse.hi[d] - 1;
                        shrunk = true;
                    }
                }
                if (inner.lo[0] > inner.hi[0] || inner.lo[1] > inner.hi[1] || inner.lo[2] > inner.hi[2]) {
                    OpmLog::warning(fmt::format("--well-refine: ring {} of zone '{}' would touch "
                                                "its parent's boundary at the grid edge and is "
                                                "dropped.", r + 1, zone.pattern));
                    continue;
                }
                if (shrunk) {
                    OpmLog::warning(fmt::format("--well-refine: ring {} of zone '{}' shrunk to "
                                                "stay inside its parent at the grid edge.",
                                                r + 1, zone.pattern));
                }
                req.parentGridName = parent->name;
                std::array<int,3> units{};
                for (int d = 0; d < 3; ++d) {
                    req.startIJK[d] = (inner.lo[d] - parent->coarse.lo[d]) * parent->unitsPerCoarse[d];
                    req.endIJK[d] = (inner.hi[d] + 1 - parent->coarse.lo[d]) * parent->unitsPerCoarse[d];
                    units[d] = parent->unitsPerCoarse[d] * zone.factor[d];
                }
                current.push_back({ inner, req.name, units });
                requests.push_back(std::move(req));
            }
            previous = std::move(current);
        }
    }
    return requests;
}

} // namespace Opm

#endif // OPM_WELL_ZONE_LGR_HPP
