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

#ifndef OPM_FLUX_BOUNDARY_HPP
#define OPM_FLUX_BOUNDARY_HPP

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/io/eclipse/FluxFile.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace Opm {

class FluxBoundary
{
public:
    struct Face {
        int interiorActiveCell = -1;
        int interiorLocalCell = -1;
        int interiorGlobalCell = -1;
        int exteriorGlobalCell = -1;
        FaceDir::DirEnum direction = FaceDir::Unknown;
        bool isNnc = false;
        double transmissibility = 0.0;

        //! \brief PVT region of the cell on the far side, as recorded by the
        //!        producing run.
        int exteriorPvtRegion = 0;

        //! \brief Depth of the centre of the cell on the far side.
        //!
        //! \details NaN if the producing run did not record it, in which case
        //!   the recorded pressure has to be taken at face value.
        double exteriorDepth = std::numeric_limits<double>::quiet_NaN();

        //! \brief Equilibration region of the cell on the far side, zero-based.
        //!
        //! \details Needed to find the threshold pressure across the face:
        //!   the sector's own EQLNUM stops at its edge, so it cannot tell what
        //!   region lies beyond. Negative if the producing run did not record
        //!   it, in which case no threshold can be applied.
        int exteriorEquilRegion = -1;

        bool operator==(const Face& other) const = default;
    };

    static std::vector<int> buildLocalToActive(const std::vector<int>& localToGlobal);

    static FluxBoundary fromData(const EclIO::FluxFile::Data& data,
                                 const std::vector<int>& localToActive);

    static FluxBoundary load(const std::string& filename,
                             const std::vector<int>& localToActive);

    static std::filesystem::path selectInputPath(const std::filesystem::path& inputDir,
                                                 const std::string& baseName);

    static const EclIO::FluxFile::ReportStep* selectReportStep(const EclIO::FluxFile::Data& data,
                                                                int episodeIndex);

    /// Select the boundary record covering \p time, in seconds.
    ///
    /// Records cover consecutive intervals \c [startTime, startTime+stepLength],
    /// so this works whether the producer wrote one record per time step or one
    /// per report step. Queries outside the recorded range clamp to the first or
    /// the last record.
    ///
    /// Returns nullptr only when there are no records at all.
    static const EclIO::FluxFile::ReportStep* selectRecordAt(const EclIO::FluxFile::Data& data,
                                                              double time);

    std::array<std::vector<int>, 6> buildDirectionalFaceIndices(std::size_t numActiveCells) const;

    /// Boundary faces that reach their neighbour through a non-neighbour
    /// connection, indexed by the interior cell they belong to.
    ///
    /// These cannot go in the directional index: they have no face direction,
    /// because the cell they connect to is not an axis neighbour, and one cell
    /// may carry several of them. The reduced run has no such connection in its
    /// own grid either, so what crosses one enters as a source in the interior
    /// cell rather than through a face.
    std::vector<std::vector<int>> buildNncFaceIndices(std::size_t numActiveCells) const;

    const Face* faceFromSlot(int slot) const;

    //! \brief Replace the default outer-boundary transmissibility of every
    //!        sector face with the one the parent run had across it.
    //!
    //! \param[in] boundaryFaceOrdinal Where a face sits in its cell's list of
    //!   boundary faces, given the interior active cell and the direction. That
    //!   position, not the direction, is how the discretisation addresses a
    //!   boundary face: the stencil appends one per intersection without a
    //!   neighbour, in iteration order, so a cell whose I- neighbour is inside
    //!   the sector carries its I+ face at position zero. Faces it reports as
    //!   absent have no grid face here and are left alone.
    template <class TransmissibilityLike, class BoundaryFaceOrdinal>
    std::size_t applyTransmissibilityOverrides(TransmissibilityLike& transmissibility,
                                               BoundaryFaceOrdinal&& boundaryFaceOrdinal) const
    {
        std::size_t appliedCount = 0;
        for (const auto& face : this->faces_) {
            if (face.isNnc || face.direction == FaceDir::Unknown) {
                continue;
            }
            if (face.interiorActiveCell < 0) {
                continue;
            }
            if (!std::isfinite(face.transmissibility)) {
                continue;
            }

            const auto ordinal = boundaryFaceOrdinal(face.interiorActiveCell, face.direction);
            if (ordinal < 0) {
                continue;
            }

            // A face the parent could not flow through must not become an open
            // boundary here. Zero is applied rather than skipped, because
            // skipping would leave the default outer-boundary value in place
            // and turn a sealing fault, or a MULTREGT of zero, into a leak.
            const auto value = (face.transmissibility > 0.0)
                ? face.transmissibility
                : decltype(face.transmissibility){0};

            transmissibility.setTransmissibilityBoundary(static_cast<unsigned>(face.interiorActiveCell),
                                                         static_cast<unsigned>(ordinal),
                                                         value);
            ++appliedCount;
        }

        return appliedCount;
    }

    const EclIO::FluxFile::Data& data() const;
    const std::vector<Face>& faces() const;

private:
    EclIO::FluxFile::Data data_;
    std::vector<Face> faces_;
};

} // namespace Opm

#endif // OPM_FLUX_BOUNDARY_HPP
