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

#include <config.h>

#define BOOST_TEST_MODULE OPM_test_FluxRegions

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/simulators/flow/flux/FluxRegions.hpp>

#include <array>
#include <stdexcept>
#include <vector>

namespace {

int globalIndex(const std::array<int, 3>& dims, const int i, const int j, const int k)
{
    return Opm::FluxRegions::cartesianIndex(dims, i - 1, j - 1, k - 1);
}

} // namespace

BOOST_AUTO_TEST_CASE(ExtractsMinimalBoxAndLocalMap)
{
    const std::array<int, 3> dims{4, 3, 2};
    std::vector<int> regions(dims[0] * dims[1] * dims[2], 0);

    regions[globalIndex(dims, 2, 2, 1)] = 2;
    regions[globalIndex(dims, 3, 2, 1)] = 2;
    regions[globalIndex(dims, 2, 3, 2)] = 2;

    const auto extracted = Opm::FluxRegions::extract(dims, regions);
    BOOST_REQUIRE_EQUAL(extracted.size(), 1U);

    const auto& region = extracted.front();
    const auto expectedBox = Opm::FluxRegions::Box{2, 3, 2, 3, 1, 2};
    BOOST_CHECK_EQUAL(region.regionId, 2);
    BOOST_CHECK(region.box == expectedBox);
    BOOST_CHECK_EQUAL(region.localToGlobal.size(), 8U);
    BOOST_CHECK_EQUAL(region.selectedGlobalCells.size(), 3U);

    BOOST_CHECK_EQUAL(region.localToGlobal[0], globalIndex(dims, 2, 2, 1));
    BOOST_CHECK_EQUAL(region.localToGlobal[1], globalIndex(dims, 3, 2, 1));
    BOOST_CHECK_EQUAL(region.localToGlobal[6], globalIndex(dims, 2, 3, 2));
    BOOST_CHECK_EQUAL(region.localToGlobal[2], -1);
    BOOST_CHECK_EQUAL(region.localToGlobal[7], -1);
}

BOOST_AUTO_TEST_CASE(ExtractsMultipleRegionsInAscendingIdOrder)
{
    const std::array<int, 3> dims{3, 2, 1};
    std::vector<int> regions(dims[0] * dims[1] * dims[2], 0);

    regions[globalIndex(dims, 1, 1, 1)] = 5;
    regions[globalIndex(dims, 3, 2, 1)] = 2;

    const auto extracted = Opm::FluxRegions::extract(dims, regions);
    BOOST_REQUIRE_EQUAL(extracted.size(), 2U);
    BOOST_CHECK_EQUAL(extracted[0].regionId, 2);
    BOOST_CHECK_EQUAL(extracted[1].regionId, 5);
}

BOOST_AUTO_TEST_CASE(RejectsWrongRegionVectorSize)
{
    const std::array<int, 3> dims{2, 2, 1};
    BOOST_CHECK_THROW(Opm::FluxRegions::extract(dims, std::vector<int>{1, 2, 3}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(DetectsUniqueSelectedRegion)
{
    const std::vector<int> regions{0, 4, 4, 0, 4};
    BOOST_CHECK_EQUAL(Opm::FluxRegions::uniqueSelectedRegion(regions), 4);
}

BOOST_AUTO_TEST_CASE(RejectsMissingOrAmbiguousSelectedRegion)
{
    BOOST_CHECK_THROW(Opm::FluxRegions::uniqueSelectedRegion(std::vector<int>{0, 0, 0}), std::invalid_argument);
    BOOST_CHECK_THROW(Opm::FluxRegions::uniqueSelectedRegion(std::vector<int>{1, 0, 2}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(BuildsActnumForSelectedRegion)
{
    const auto actnum = Opm::FluxRegions::buildActnum(std::vector<int>{0, 2, 3, 2, 0}, 2);
    BOOST_REQUIRE_EQUAL(actnum.size(), 5U);
    BOOST_CHECK_EQUAL(actnum[0], 0);
    BOOST_CHECK_EQUAL(actnum[1], 1);
    BOOST_CHECK_EQUAL(actnum[2], 0);
    BOOST_CHECK_EQUAL(actnum[3], 1);
    BOOST_CHECK_EQUAL(actnum[4], 0);
}

BOOST_AUTO_TEST_CASE(BuildsActnumForSelectedRegionList)
{
    const auto actnum = Opm::FluxRegions::buildActnum(std::vector<int>{0, 2, 3, 2, 4}, std::vector<int>{2, 4});
    BOOST_REQUIRE_EQUAL(actnum.size(), 5U);
    BOOST_CHECK_EQUAL(actnum[0], 0);
    BOOST_CHECK_EQUAL(actnum[1], 1);
    BOOST_CHECK_EQUAL(actnum[2], 0);
    BOOST_CHECK_EQUAL(actnum[3], 1);
    BOOST_CHECK_EQUAL(actnum[4], 1);
}

BOOST_AUTO_TEST_CASE(ExtractsCartesianBoundaryFaces)
{
    const std::array<int, 3> dims{4, 1, 1};
    std::vector<int> regions(dims[0] * dims[1] * dims[2], 0);

    regions[globalIndex(dims, 2, 1, 1)] = 9;
    regions[globalIndex(dims, 3, 1, 1)] = 9;

    const auto extracted = Opm::FluxRegions::extract(dims, regions);
    BOOST_REQUIRE_EQUAL(extracted.size(), 1U);

    const auto& region = extracted.front();
    BOOST_REQUIRE_EQUAL(region.boundaryFaces.size(), 2U);

    const auto left = Opm::FluxRegions::BoundaryFace{
        0,
        globalIndex(dims, 2, 1, 1),
        globalIndex(dims, 1, 1, 1),
        Opm::FaceDir::XMinus,
        false,
    };
    const auto right = Opm::FluxRegions::BoundaryFace{
        1,
        globalIndex(dims, 3, 1, 1),
        globalIndex(dims, 4, 1, 1),
        Opm::FaceDir::XPlus,
        false,
    };

    BOOST_CHECK(region.boundaryFaces[0] == left);
    BOOST_CHECK(region.boundaryFaces[1] == right);
}

BOOST_AUTO_TEST_CASE(ExtractsNncBoundaryFaces)
{
    const std::array<int, 3> dims{3, 1, 1};
    std::vector<int> regions(dims[0] * dims[1] * dims[2], 0);

    regions[globalIndex(dims, 1, 1, 1)] = 7;

    const std::vector<std::array<int, 2>> nnc{
        {globalIndex(dims, 1, 1, 1), globalIndex(dims, 3, 1, 1)}
    };

    const auto extracted = Opm::FluxRegions::extract(dims, regions, nnc);
    BOOST_REQUIRE_EQUAL(extracted.size(), 1U);

    const auto& region = extracted.front();
    BOOST_REQUIRE_EQUAL(region.boundaryFaces.size(), 2U);

    const auto cartesianFace = Opm::FluxRegions::BoundaryFace{
        0,
        globalIndex(dims, 1, 1, 1),
        globalIndex(dims, 2, 1, 1),
        Opm::FaceDir::XPlus,
        false,
    };
    const auto nncFace = Opm::FluxRegions::BoundaryFace{
        0,
        globalIndex(dims, 1, 1, 1),
        globalIndex(dims, 3, 1, 1),
        Opm::FaceDir::Unknown,
        true,
    };

    BOOST_CHECK(region.boundaryFaces[0] == cartesianFace);
    BOOST_CHECK(region.boundaryFaces[1] == nncFace);
}

// The region map covers every cell of the grid, including the inactive ones.
// An inactive cell holds no fluid and takes no part in the flow, so it must not
// become part of the region: the consumer maps boundary faces onto active cells
// by counting the entries of localToGlobal, and counting a cell that the sector
// grid does not have would shift every later face onto the wrong cell.
BOOST_AUTO_TEST_CASE(InactiveCellsAreExcludedFromTheRegion)
{
    const std::array<int, 3> dims{4, 1, 1};
    std::vector<int> regions(dims[0] * dims[1] * dims[2], 0);
    regions[globalIndex(dims, 2, 1, 1)] = 1;
    regions[globalIndex(dims, 3, 1, 1)] = 1;

    // Without ACTNUM both cells belong to the region.
    {
        const auto extracted = Opm::FluxRegions::extract(dims, regions);
        BOOST_REQUIRE_EQUAL(extracted.size(), 1U);
        BOOST_CHECK_EQUAL(extracted.front().selectedGlobalCells.size(), 2U);
        BOOST_CHECK_EQUAL(extracted.front().localToGlobal.size(), 2U);
    }

    // Marking the second of them inactive leaves a one-cell region.
    std::vector<int> actnum(dims[0] * dims[1] * dims[2], 1);
    actnum[globalIndex(dims, 3, 1, 1)] = 0;

    const auto extracted = Opm::FluxRegions::extract(dims, regions, actnum, {});
    BOOST_REQUIRE_EQUAL(extracted.size(), 1U);

    const auto& region = extracted.front();
    BOOST_REQUIRE_EQUAL(region.selectedGlobalCells.size(), 1U);
    BOOST_CHECK_EQUAL(region.selectedGlobalCells.front(), globalIndex(dims, 2, 1, 1));

    // The bounding box shrinks with the region, so localToGlobal holds exactly
    // the one active cell and nothing else.
    const auto expectedBox = Opm::FluxRegions::Box{2, 2, 1, 1, 1, 1};
    BOOST_CHECK(region.box == expectedBox);
    BOOST_REQUIRE_EQUAL(region.localToGlobal.size(), 1U);
    BOOST_CHECK_EQUAL(region.localToGlobal.front(), globalIndex(dims, 2, 1, 1));
}

// Nothing flows across a face towards an inactive cell, so that face is a
// no-flow boundary rather than one the parent run has to supply data for.
BOOST_AUTO_TEST_CASE(FacesTowardsInactiveCellsAreNotBoundaryFaces)
{
    const std::array<int, 3> dims{3, 1, 1};
    std::vector<int> regions(dims[0] * dims[1] * dims[2], 0);
    regions[globalIndex(dims, 2, 1, 1)] = 1;

    // With every cell active the single region cell has two boundary faces.
    {
        const auto extracted = Opm::FluxRegions::extract(dims, regions);
        BOOST_REQUIRE_EQUAL(extracted.size(), 1U);
        BOOST_CHECK_EQUAL(extracted.front().boundaryFaces.size(), 2U);
    }

    // Deactivating the neighbour on one side removes that face.
    std::vector<int> actnum(dims[0] * dims[1] * dims[2], 1);
    actnum[globalIndex(dims, 1, 1, 1)] = 0;

    const auto extracted = Opm::FluxRegions::extract(dims, regions, actnum, {});
    BOOST_REQUIRE_EQUAL(extracted.size(), 1U);

    const auto& faces = extracted.front().boundaryFaces;
    BOOST_REQUIRE_EQUAL(faces.size(), 1U);
    BOOST_CHECK_EQUAL(faces.front().exteriorGlobalCell, globalIndex(dims, 3, 1, 1));
}

BOOST_AUTO_TEST_CASE(RejectsActnumOfWrongSize)
{
    const std::array<int, 3> dims{3, 1, 1};
    std::vector<int> regions(dims[0] * dims[1] * dims[2], 0);
    regions[globalIndex(dims, 2, 1, 1)] = 1;

    const std::vector<int> actnum(2, 1);
    BOOST_CHECK_THROW(Opm::FluxRegions::extract(dims, regions, actnum, {}),
                      std::invalid_argument);
}