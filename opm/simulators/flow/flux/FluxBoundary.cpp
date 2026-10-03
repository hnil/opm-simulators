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

#include <opm/simulators/flow/flux/FluxBoundary.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <stdexcept>

namespace Opm {

//! \details Numbers the region's cells in the order they appear, which is the
//!   order a serial reduced run activates them in, so the two agree. A
//!   parallel run does not: each rank holds an arbitrary subset of the region
//!   under its own numbering, and its caller has to derive the mapping from
//!   the grid instead of from the file.
std::vector<int> FluxBoundary::buildLocalToActive(const std::vector<int>& localToGlobal)
{
    std::vector<int> localToActive(localToGlobal.size(), -1);
    int activeIndex = 0;
    for (std::size_t localIndex = 0; localIndex < localToGlobal.size(); ++localIndex) {
        if (localToGlobal[localIndex] >= 0) {
            localToActive[localIndex] = activeIndex;
            ++activeIndex;
        }
    }

    return localToActive;
}

FluxBoundary FluxBoundary::fromData(const EclIO::FluxFile::Data& data,
                                    const std::vector<int>& localToActive)
{
    if (localToActive.size() != data.localToGlobal.size()) {
        OPM_THROW(std::invalid_argument,
                  "FluxBoundary: localToActive size must match FLUX localToGlobal size");
    }

    FluxBoundary boundary;
    boundary.data_ = data;
    boundary.faces_.reserve(data.boundaryFaces.size());

    for (const auto& face : data.boundaryFaces) {
        if (face.interiorLocalCell < 0
            || static_cast<std::size_t>(face.interiorLocalCell) >= localToActive.size()) {
            OPM_THROW(std::invalid_argument,
                      "FluxBoundary: boundary face references invalid interior local cell");
        }

        // A face whose interior cell this rank does not hold keeps its place in
        // the list with no active cell of its own. Dropping it would renumber
        // everything after it, and the file's payload arrays are indexed by
        // face position, so every rank carries the whole list whether or not it
        // can act on a given face.
        const auto interiorActive = localToActive[face.interiorLocalCell];

        const auto interiorGlobal = data.localToGlobal[face.interiorLocalCell];
        if (interiorGlobal < 0) {
            OPM_THROW(std::invalid_argument,
                      "FluxBoundary: localToGlobal missing interior global cell for boundary face");
        }

        boundary.faces_.push_back(Face{
            interiorActive,
            face.interiorLocalCell,
            interiorGlobal,
            face.exteriorGlobalCell,
            static_cast<FaceDir::DirEnum>(face.direction),
            face.direction == static_cast<int>(FaceDir::Unknown),
            face.transmissibility,
            face.exteriorPvtRegion,
            face.exteriorDepth,
            face.exteriorEquilRegion,
        });
    }

    return boundary;
}

FluxBoundary FluxBoundary::load(const std::string& filename,
                                const std::vector<int>& localToActive)
{
    return fromData(EclIO::FluxFile::read(filename), localToActive);
}

std::filesystem::path FluxBoundary::selectInputPath(const std::filesystem::path& inputDir,
                                                    const std::string& baseName)
{
    const auto fluxPath = inputDir / (baseName + ".FLUX");
    if (std::filesystem::exists(fluxPath)) {
        return fluxPath;
    }

    const auto prefix = baseName + ".FLUX";
    std::vector<std::filesystem::path> suffixCandidates;
    if (std::filesystem::exists(inputDir)) {
        for (const auto& entry : std::filesystem::directory_iterator(inputDir)) {
            if (!entry.is_regular_file()) {
                continue;
            }

            const auto filename = entry.path().filename().string();
            if (filename.rfind(prefix, 0) != 0) {
                continue;
            }

            const auto suffix = filename.substr(prefix.size());
            if (suffix.size() != 4
                || !std::all_of(suffix.begin(), suffix.end(), [](const char c)
                                {
                                    return std::isdigit(static_cast<unsigned char>(c)) != 0;
                                })) {
                continue;
            }

            suffixCandidates.push_back(entry.path());
        }
    }

    if (!suffixCandidates.empty()) {
        std::sort(suffixCandidates.begin(), suffixCandidates.end());
        return suffixCandidates.front();
    }

    return fluxPath;
}

const EclIO::FluxFile::ReportStep* FluxBoundary::selectReportStep(const EclIO::FluxFile::Data& data,
                                                                   const int episodeIndex)
{
    if (data.reportSteps.empty()) {
        return nullptr;
    }

    // The producer appends a report step using reportStepNum = episodeIndex + 1,
    // and the stored rates are averaged over that episode. A consumer running
    // episode e therefore has to use the entry with reportStep == e + 1.
    const auto byOneBased = episodeIndex + 1;
    const auto byOneBasedCount = static_cast<int>(std::count_if(data.reportSteps.begin(),
                                                                 data.reportSteps.end(),
                                                                 [byOneBased](const auto& step)
                                                                 {
                                                                     return step.reportStep == byOneBased;
                                                                 }));
    if (byOneBasedCount > 0) {
        const auto it = std::find_if(data.reportSteps.begin(),
                                     data.reportSteps.end(),
                                     [byOneBased](const auto& step)
                                     {
                                         return step.reportStep == byOneBased;
                                     });
        return &(*it);
    }

    const auto byReportStep = static_cast<int>(std::count_if(data.reportSteps.begin(),
                                                              data.reportSteps.end(),
                                                              [episodeIndex](const auto& step)
                                                              {
                                                                  return step.reportStep == episodeIndex;
                                                              }));
    if (byReportStep > 0) {
        const auto it = std::find_if(data.reportSteps.begin(),
                                     data.reportSteps.end(),
                                     [episodeIndex](const auto& step)
                                     {
                                         return step.reportStep == episodeIndex;
                                     });
        return &(*it);
    }

    const auto nonNegativeEpisode = std::max(episodeIndex, 0);
    const auto index = std::min(static_cast<std::size_t>(nonNegativeEpisode),
                                data.reportSteps.size() - 1);
    return &data.reportSteps[index];
}

const EclIO::FluxFile::ReportStep* FluxBoundary::selectRecordAt(const EclIO::FluxFile::Data& data,
                                                                 const double time)
{
    if (data.reportSteps.empty()) {
        return nullptr;
    }

    // A record covers (startTime, startTime + stepLength], and 'time' is the
    // start of the step that is about to be taken. The record that applies is
    // therefore the first one that ends strictly after it: a record ending
    // exactly at 'time' describes flow that has already happened.
    const auto pos = std::lower_bound(data.reportSteps.begin(), data.reportSteps.end(), time,
                                      [](const auto& step, const double t)
                                      {
                                          return (step.startTime + step.stepLength) <= t;
                                      });

    if (pos == data.reportSteps.end()) {
        // After the last record.
        return &data.reportSteps.back();
    }

    return &(*pos);
}

std::array<std::vector<int>, 6> FluxBoundary::buildDirectionalFaceIndices(const std::size_t numActiveCells) const
{
    std::array<std::vector<int>, 6> directionalFaceIndices;
    for (auto& values : directionalFaceIndices) {
        values.assign(numActiveCells, 0);
    }

    const auto dirToIndex = [](const FaceDir::DirEnum dir)
    {
        int idx = 0;
        int div = static_cast<int>(dir);
        while ((div /= 2) >= 1) {
            ++idx;
        }
        return idx;
    };

    for (std::size_t faceIndex = 0; faceIndex < this->faces_.size(); ++faceIndex) {
        const auto& face = this->faces_[faceIndex];
        if (face.isNnc || face.direction == FaceDir::Unknown) {
            continue;
        }
        // Not this rank's face to impose. In parallel most of the region's
        // boundary belongs to somebody else.
        if (face.interiorActiveCell < 0) {
            continue;
        }
        if (static_cast<std::size_t>(face.interiorActiveCell) >= numActiveCells) {
            OPM_THROW(std::invalid_argument,
                      "FluxBoundary: face references invalid active cell for directional mapping");
        }

        auto& slot = directionalFaceIndices[dirToIndex(face.direction)][face.interiorActiveCell];
        if (slot != 0) {
            OPM_THROW(std::invalid_argument,
                      "FluxBoundary: duplicate FLUX boundary face for active cell and direction");
        }

        slot = static_cast<int>(faceIndex) + 1;
    }

    return directionalFaceIndices;
}

std::vector<std::vector<int>>
FluxBoundary::buildNncFaceIndices(const std::size_t numActiveCells) const
{
    std::vector<std::vector<int>> nncFaceIndices(numActiveCells);

    for (std::size_t faceIndex = 0; faceIndex < this->faces_.size(); ++faceIndex) {
        const auto& face = this->faces_[faceIndex];
        if (!face.isNnc && (face.direction != FaceDir::Unknown)) {
            continue;
        }

        if (face.interiorActiveCell < 0) {
            continue;
        }
        if (static_cast<std::size_t>(face.interiorActiveCell) >= numActiveCells) {
            OPM_THROW(std::invalid_argument,
                      "FluxBoundary: face references invalid active cell for NNC mapping");
        }

        // Unlike the directional faces these are not unique per cell: a cell
        // can sit at the end of several non-neighbour connections.
        nncFaceIndices[face.interiorActiveCell].push_back(static_cast<int>(faceIndex));
    }

    return nncFaceIndices;
}

const FluxBoundary::Face* FluxBoundary::faceFromSlot(const int slot) const
{
    if (slot <= 0) {
        return nullptr;
    }

    const auto faceIndex = static_cast<std::size_t>(slot - 1);
    if (faceIndex >= this->faces_.size()) {
        return nullptr;
    }

    return &this->faces_[faceIndex];
}

const EclIO::FluxFile::Data& FluxBoundary::data() const
{
    return this->data_;
}

const std::vector<FluxBoundary::Face>& FluxBoundary::faces() const
{
    return this->faces_;
}

} // namespace Opm
