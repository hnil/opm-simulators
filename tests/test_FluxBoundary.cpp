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

#define BOOST_TEST_MODULE OPM_test_FluxBoundary

#include <boost/test/unit_test.hpp>

#include <opm/io/eclipse/FluxFile.hpp>
#include <opm/simulators/flow/flux/FluxBoundary.hpp>
#include <opm/simulators/flow/flux/FluxDumper.hpp>
#include <opm/simulators/flow/flux/FluxRegions.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace {

int globalIndex(const std::array<int, 3>& dims, const int i, const int j, const int k)
{
    return Opm::FluxRegions::cartesianIndex(dims, i - 1, j - 1, k - 1);
}

} // namespace

BOOST_AUTO_TEST_CASE(BuildsContiguousLocalToActiveMap)
{
    const std::vector<int> localToGlobal{10, -1, 12, 17, -1};
    const auto localToActive = Opm::FluxBoundary::buildLocalToActive(localToGlobal);

    BOOST_REQUIRE_EQUAL(localToActive.size(), 5U);
    BOOST_CHECK_EQUAL(localToActive[0], 0);
    BOOST_CHECK_EQUAL(localToActive[1], -1);
    BOOST_CHECK_EQUAL(localToActive[2], 1);
    BOOST_CHECK_EQUAL(localToActive[3], 2);
    BOOST_CHECK_EQUAL(localToActive[4], -1);
}

BOOST_AUTO_TEST_CASE(MapsBoundaryFacesFromFluxFileData)
{
    const std::array<int, 3> dims{4, 1, 1};
    std::vector<int> regionValues(dims[0] * dims[1] * dims[2], 0);
    regionValues[globalIndex(dims, 2, 1, 1)] = 3;
    regionValues[globalIndex(dims, 3, 1, 1)] = 3;

    const auto regions = Opm::FluxRegions::extract(dims, regionValues);
    BOOST_REQUIRE_EQUAL(regions.size(), 1U);

    Opm::FluxDumper dumper("PARENT", 3, dims, regions.front(),
                           Opm::EclIO::FluxFile::Mode::Flux,
                           Opm::EclIO::FluxFile::Sampling::Instant,
                           static_cast<int>(Opm::EclIO::FluxFile::Phase::Oil));

    Opm::FluxDumper::ReportStepData step;
    step.reportStep = 0;
    step.simStep = 0;
    step.startTime = 0.0;
    step.stepLength = 1.0;
    step.massRates = {1.0, 2.0};
    dumper.appendReportStep(step);

    const auto localToActive = Opm::FluxBoundary::buildLocalToActive(regions.front().localToGlobal);
    const auto boundary = Opm::FluxBoundary::fromData(dumper.data(), localToActive);

    BOOST_REQUIRE_EQUAL(boundary.faces().size(), 2U);
    BOOST_CHECK_EQUAL(boundary.faces()[0].interiorActiveCell, 0);
    BOOST_CHECK_EQUAL(boundary.faces()[0].interiorLocalCell, 0);
    BOOST_CHECK_EQUAL(boundary.faces()[0].interiorGlobalCell, globalIndex(dims, 2, 1, 1));
    BOOST_CHECK_EQUAL(boundary.faces()[0].direction, Opm::FaceDir::XMinus);
    BOOST_CHECK(!boundary.faces()[0].isNnc);

    BOOST_CHECK_EQUAL(boundary.faces()[1].interiorActiveCell, 1);
    BOOST_CHECK_EQUAL(boundary.faces()[1].interiorLocalCell, 1);
    BOOST_CHECK_EQUAL(boundary.faces()[1].interiorGlobalCell, globalIndex(dims, 3, 1, 1));
    BOOST_CHECK_EQUAL(boundary.faces()[1].direction, Opm::FaceDir::XPlus);
}

BOOST_AUTO_TEST_CASE(LoadsFluxFileAndPreservesNncFaces)
{
    const std::array<int, 3> dims{3, 1, 1};
    std::vector<int> regionValues(dims[0] * dims[1] * dims[2], 0);
    regionValues[globalIndex(dims, 1, 1, 1)] = 8;

    const std::vector<std::array<int, 2>> nnc{{globalIndex(dims, 1, 1, 1), globalIndex(dims, 3, 1, 1)}};
    const auto regions = Opm::FluxRegions::extract(dims, regionValues, nnc);
    BOOST_REQUIRE_EQUAL(regions.size(), 1U);

    Opm::FluxDumper dumper("PARENT", 8, dims, regions.front(),
                           Opm::EclIO::FluxFile::Mode::Flux,
                           Opm::EclIO::FluxFile::Sampling::Instant,
                           static_cast<int>(Opm::EclIO::FluxFile::Phase::Oil));

    Opm::FluxDumper::ReportStepData step;
    step.reportStep = 0;
    step.simStep = 0;
    step.startTime = 0.0;
    step.stepLength = 1.0;
    step.massRates = {11.0, 99.0};
    dumper.appendReportStep(step);

    const auto outPath = std::filesystem::path{"test_fluxboundary_roundtrip.FLUX"};
    dumper.flush(outPath.string(), false);

    const auto localToActive = Opm::FluxBoundary::buildLocalToActive(regions.front().localToGlobal);
    const auto boundary = Opm::FluxBoundary::load(outPath.string(), localToActive);

    BOOST_REQUIRE_EQUAL(boundary.faces().size(), 2U);
    BOOST_CHECK(!boundary.faces()[0].isNnc);
    BOOST_CHECK(boundary.faces()[1].isNnc);
    BOOST_CHECK_EQUAL(boundary.faces()[1].direction, Opm::FaceDir::Unknown);
    BOOST_CHECK_EQUAL(boundary.faces()[1].exteriorGlobalCell, globalIndex(dims, 3, 1, 1));

    std::filesystem::remove(outPath);
}

BOOST_AUTO_TEST_CASE(RejectsMismatchedLocalToActiveSize)
{
    Opm::EclIO::FluxFile::Data data;
    data.localToGlobal = {10};
    data.boundaryFaces = {{0, static_cast<int>(Opm::FaceDir::XMinus), 9, 0.0}};

    BOOST_CHECK_THROW(Opm::FluxBoundary::fromData(data, {}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(KeepsBoundaryFaceWhoseInteriorCellIsElsewhere)
{
    // In parallel a rank holds only part of the region, so most of the
    // boundary belongs to somebody else. Such a face keeps its place in the
    // list, because the file's payload arrays are indexed by face position and
    // dropping one would renumber every face after it. It carries no active
    // cell, and the lookups a rank uses to impose a boundary pass over it.
    Opm::EclIO::FluxFile::Data data;
    data.localToGlobal = {10, 11};
    data.boundaryFaces = {
        {0, static_cast<int>(Opm::FaceDir::XMinus), 9, 1.5},
        {1, static_cast<int>(Opm::FaceDir::XPlus), 12, 2.5},
    };

    // The first cell is on another rank; only the second is here.
    const auto boundary = Opm::FluxBoundary::fromData(data, {-1, 0});

    BOOST_REQUIRE_EQUAL(boundary.faces().size(), 2U);
    BOOST_CHECK_EQUAL(boundary.faces()[0].interiorActiveCell, -1);
    BOOST_CHECK_EQUAL(boundary.faces()[1].interiorActiveCell, 0);

    // Face positions are unchanged, so the second face still resolves to slot
    // 2 and indexes the second entry of every payload array. The direction
    // arrays are keyed as elsewhere in the flux code, by the bit position of
    // the direction: 0 is I+ and 1 is I-.
    const auto directional = boundary.buildDirectionalFaceIndices(1);
    BOOST_CHECK_EQUAL(directional[0][0], 2);
    BOOST_CHECK_EQUAL(directional[1][0], 0);
}

BOOST_AUTO_TEST_CASE(BuildsDirectionalFaceIndexRegistration)
{
    const std::array<int, 3> dims{4, 1, 1};
    std::vector<int> regionValues(dims[0] * dims[1] * dims[2], 0);
    regionValues[globalIndex(dims, 2, 1, 1)] = 3;
    regionValues[globalIndex(dims, 3, 1, 1)] = 3;

    const auto regions = Opm::FluxRegions::extract(dims, regionValues);
    BOOST_REQUIRE_EQUAL(regions.size(), 1U);

    Opm::FluxDumper dumper("PARENT", 3, dims, regions.front(),
                           Opm::EclIO::FluxFile::Mode::Flux,
                           Opm::EclIO::FluxFile::Sampling::Instant,
                           static_cast<int>(Opm::EclIO::FluxFile::Phase::Oil));

    Opm::FluxDumper::ReportStepData step;
    step.reportStep = 0;
    step.simStep = 0;
    step.startTime = 0.0;
    step.stepLength = 1.0;
    step.massRates = {1.0, 2.0};
    dumper.appendReportStep(step);

    const auto localToActive = Opm::FluxBoundary::buildLocalToActive(regions.front().localToGlobal);
    const auto boundary = Opm::FluxBoundary::fromData(dumper.data(), localToActive);
    const auto directional = boundary.buildDirectionalFaceIndices(2);

    BOOST_CHECK_EQUAL(directional[1][0], 1);
    BOOST_CHECK_EQUAL(directional[0][1], 2);
    BOOST_CHECK_EQUAL(directional[0][0], 0);
    BOOST_CHECK_EQUAL(directional[1][1], 0);
}

BOOST_AUTO_TEST_CASE(RejectsDuplicateDirectionalFaceRegistration)
{
    Opm::EclIO::FluxFile::Data data;
    data.localToGlobal = {10};
    data.boundaryFaces = {
        {0, static_cast<int>(Opm::FaceDir::XMinus), 9, 0.0},
        {0, static_cast<int>(Opm::FaceDir::XMinus), 8, 0.0},
    };

    const auto boundary = Opm::FluxBoundary::fromData(data, {0});
    BOOST_CHECK_THROW(boundary.buildDirectionalFaceIndices(1), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(ReturnsFaceFromDirectionalSlot)
{
    Opm::EclIO::FluxFile::Data data;
    data.localToGlobal = {10, 11};
    data.boundaryFaces = {
        {0, static_cast<int>(Opm::FaceDir::XMinus), 9, 1.0},
        {1, static_cast<int>(Opm::FaceDir::XPlus), 12, 2.0},
    };

    const auto boundary = Opm::FluxBoundary::fromData(data, {0, 1});
    BOOST_REQUIRE(boundary.faceFromSlot(1) != nullptr);
    BOOST_REQUIRE(boundary.faceFromSlot(2) != nullptr);
    BOOST_CHECK_EQUAL(boundary.faceFromSlot(1)->interiorActiveCell, 0);
    BOOST_CHECK_EQUAL(boundary.faceFromSlot(2)->interiorActiveCell, 1);
    BOOST_CHECK(boundary.faceFromSlot(0) == nullptr);
    BOOST_CHECK(boundary.faceFromSlot(3) == nullptr);
}

BOOST_AUTO_TEST_CASE(SelectsReportStepByEpisodeThenFallback)
{
    Opm::EclIO::FluxFile::Data data;
    const auto mkStep = [](const int reportStep, const int simStep)
    {
        Opm::EclIO::FluxFile::ReportStep step;
        step.reportStep = reportStep;
        step.simStep = simStep;
        step.startTime = 0.0;
        step.stepLength = 0.1;
        return step;
    };
    data.reportSteps = {
        mkStep(2, 20),
        mkStep(4, 40),
        mkStep(8, 80),
    };

    const auto* exact = Opm::FluxBoundary::selectReportStep(data, 4);
    BOOST_REQUIRE(exact != nullptr);
    BOOST_CHECK_EQUAL(exact->reportStep, 4);

    const auto* oneBased = Opm::FluxBoundary::selectReportStep(data, 1);
    BOOST_REQUIRE(oneBased != nullptr);
    BOOST_CHECK_EQUAL(oneBased->reportStep, 2);

    const auto* clamped = Opm::FluxBoundary::selectReportStep(data, 99);
    BOOST_REQUIRE(clamped != nullptr);
    BOOST_CHECK_EQUAL(clamped->reportStep, 8);

    Opm::EclIO::FluxFile::Data empty;
    BOOST_CHECK(Opm::FluxBoundary::selectReportStep(empty, 0) == nullptr);
}

BOOST_AUTO_TEST_CASE(AppliesTransmissibilityOverridesWithFiltering)
{
    struct DummyTransmissibility {
        std::vector<std::tuple<unsigned, unsigned, double>> calls;

        void setTransmissibilityBoundary(const unsigned elemIdx,
                                         const unsigned boundaryFaceIdx,
                                         const double value)
        {
            this->calls.emplace_back(elemIdx, boundaryFaceIdx, value);
        }
    };

    Opm::EclIO::FluxFile::Data data;
    data.localToGlobal = {10, 11, 12, 13, 14, 15};
    data.boundaryFaces = {
        {0, static_cast<int>(Opm::FaceDir::XMinus), 9, 4.5},
        {1, static_cast<int>(Opm::FaceDir::XPlus), 12, 7.25},
        {2, static_cast<int>(Opm::FaceDir::Unknown), 99, 3.0},
        {3, static_cast<int>(Opm::FaceDir::YMinus), 8, 0.0},
        {4, static_cast<int>(Opm::FaceDir::YPlus), 7, -2.0},
        {5, static_cast<int>(Opm::FaceDir::ZMinus), 6, std::numeric_limits<double>::quiet_NaN()},
    };

    const auto boundary = Opm::FluxBoundary::fromData(data, {0, 1, 2, 3, 4, 5});
    DummyTransmissibility transmissibility;

    // Where a face sits in its cell's boundary-face list, which is how the
    // discretisation addresses it. Deliberately not the direction index: cell 1
    // carries its I+ face at position zero, which is what happens when the
    // cell's I- neighbour lies inside the sector. Cell 4 has no grid face at
    // all and is reported absent.
    const auto ordinal = [](const int activeCell, const Opm::FaceDir::DirEnum)
    {
        switch (activeCell) {
        case 0:  return 0;
        case 1:  return 0;
        case 3:  return 2;
        case 4:  return -1;
        default: return 5;
        }
    };

    const auto applied = boundary.applyTransmissibilityOverrides(transmissibility, ordinal);

    // Faces without a direction, and faces whose stored transmissibility is not
    // a number, carry no information and are left alone. Nor is there anything
    // to do for a face with no grid face in this run. A stored value of zero or
    // less does carry information -- the parent could not flow through the
    // face -- so it is applied as zero to keep the face closed.
    BOOST_CHECK_EQUAL(applied, 3U);
    BOOST_REQUIRE_EQUAL(transmissibility.calls.size(), 3U);

    BOOST_CHECK_EQUAL(std::get<0>(transmissibility.calls[0]), 0U);
    BOOST_CHECK_EQUAL(std::get<1>(transmissibility.calls[0]), 0U);
    BOOST_CHECK_CLOSE(std::get<2>(transmissibility.calls[0]), 4.5, 1e-12);

    // The one that matters: addressed by its position, not by XPlus.
    BOOST_CHECK_EQUAL(std::get<0>(transmissibility.calls[1]), 1U);
    BOOST_CHECK_EQUAL(std::get<1>(transmissibility.calls[1]), 0U);
    BOOST_CHECK_CLOSE(std::get<2>(transmissibility.calls[1]), 7.25, 1e-12);

    // The zero face is clamped to zero; the negative one never gets that far,
    // having no grid face to be applied to.
    BOOST_CHECK_EQUAL(std::get<0>(transmissibility.calls[2]), 3U);
    BOOST_CHECK_EQUAL(std::get<1>(transmissibility.calls[2]), 2U);
    BOOST_CHECK_EQUAL(std::get<2>(transmissibility.calls[2]), 0.0);
}

BOOST_AUTO_TEST_CASE(SelectInputPathPrefersExactFluxFile)
{
    const auto dir = std::filesystem::path{"test_fluxboundary_input_select_exact"};
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);

    const std::string base = "CASE";
    const auto exact = dir / "CASE.FLUX";
    const auto s2 = dir / "CASE.FLUX0002";
    const auto s7 = dir / "CASE.FLUX0007";

    std::ofstream(exact.string()).put('\n');
    std::ofstream(s2.string()).put('\n');
    std::ofstream(s7.string()).put('\n');

    const auto selected = Opm::FluxBoundary::selectInputPath(dir, base);
    BOOST_CHECK_EQUAL(selected, exact);

    std::filesystem::remove_all(dir);
}

BOOST_AUTO_TEST_CASE(SelectInputPathChoosesLowestNumericSuffix)
{
    const auto dir = std::filesystem::path{"test_fluxboundary_input_select_suffix"};
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);

    const std::string base = "CASE";
    const auto s2 = dir / "CASE.FLUX0002";
    const auto s7 = dir / "CASE.FLUX0007";
    const auto badAlpha = dir / "CASE.FLUXABCD";
    const auto badLen = dir / "CASE.FLUX00021";

    std::ofstream(s7.string()).put('\n');
    std::ofstream(s2.string()).put('\n');
    std::ofstream(badAlpha.string()).put('\n');
    std::ofstream(badLen.string()).put('\n');

    const auto selected = Opm::FluxBoundary::selectInputPath(dir, base);
    BOOST_CHECK_EQUAL(selected, s2);

    std::filesystem::remove_all(dir);
}

BOOST_AUTO_TEST_CASE(SelectInputPathReturnsDefaultWhenNoCandidatesExist)
{
    const auto dir = std::filesystem::path{"test_fluxboundary_input_select_none"};
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);

    const std::string base = "CASE";
    std::ofstream((dir / "OTHER.FLUX0001").string()).put('\n');
    std::ofstream((dir / "CASE.FLUXABC1").string()).put('\n');

    const auto selected = Opm::FluxBoundary::selectInputPath(dir, base);
    BOOST_CHECK_EQUAL(selected, dir / "CASE.FLUX");

    std::filesystem::remove_all(dir);
}

BOOST_AUTO_TEST_CASE(SelectInputPathIgnoresDirectoryCandidates)
{
    const auto dir = std::filesystem::path{"test_fluxboundary_input_select_dir"};
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);

    const std::string base = "CASE";
    std::filesystem::create_directories(dir / "CASE.FLUX0001");
    std::ofstream((dir / "CASE.FLUX0003").string()).put('\n');

    const auto selected = Opm::FluxBoundary::selectInputPath(dir, base);
    BOOST_CHECK_EQUAL(selected, dir / "CASE.FLUX0003");

    std::filesystem::remove_all(dir);
}

BOOST_AUTO_TEST_CASE(SelectInputPathPrefersLowestValidAmongMixedCandidates)
{
    const auto dir = std::filesystem::path{"test_fluxboundary_input_select_mixed"};
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);

    const std::string base = "CASE";
    std::ofstream((dir / "CASE.FLUX0007").string()).put('\n');
    std::ofstream((dir / "CASE.FLUX0002").string()).put('\n');
    std::ofstream((dir / "CASE.FLUXABCD").string()).put('\n');

    const auto selected = Opm::FluxBoundary::selectInputPath(dir, base);
    BOOST_CHECK_EQUAL(selected, dir / "CASE.FLUX0002");

    std::filesystem::remove_all(dir);
}

namespace {

// Three consecutive records covering (0,10], (10,30] and (30,60] seconds.
Opm::EclIO::FluxFile::Data makeRecordSeries()
{
    Opm::EclIO::FluxFile::Data data;

    data.reportSteps.resize(3);

    data.reportSteps[0].reportStep = 1;
    data.reportSteps[0].startTime = 0.0;
    data.reportSteps[0].stepLength = 10.0;

    data.reportSteps[1].reportStep = 1;
    data.reportSteps[1].startTime = 10.0;
    data.reportSteps[1].stepLength = 20.0;

    data.reportSteps[2].reportStep = 2;
    data.reportSteps[2].startTime = 30.0;
    data.reportSteps[2].stepLength = 30.0;

    return data;
}

} // namespace

BOOST_AUTO_TEST_CASE(SelectsRecordCoveringTheQueryTime)
{
    const auto data = makeRecordSeries();

    // Strictly inside a record.
    BOOST_CHECK_EQUAL(Opm::FluxBoundary::selectRecordAt(data, 5.0)->startTime, 0.0);
    BOOST_CHECK_EQUAL(Opm::FluxBoundary::selectRecordAt(data, 20.0)->startTime, 10.0);
    BOOST_CHECK_EQUAL(Opm::FluxBoundary::selectRecordAt(data, 45.0)->startTime, 30.0);
}

BOOST_AUTO_TEST_CASE(RecordBoundaryBelongsToTheRecordThatStartsThere)
{
    const auto data = makeRecordSeries();

    // A query exactly on a boundary must resolve to the record that BEGINS
    // there, not the one that ends there. Time steps are queried by their
    // start time, so getting this wrong shifts every step onto the previous
    // record and replays flow that has already happened.
    BOOST_CHECK_EQUAL(Opm::FluxBoundary::selectRecordAt(data, 10.0)->startTime, 10.0);
    BOOST_CHECK_EQUAL(Opm::FluxBoundary::selectRecordAt(data, 30.0)->startTime, 30.0);
    BOOST_CHECK_EQUAL(Opm::FluxBoundary::selectRecordAt(data, 60.0)->startTime, 30.0);
}

BOOST_AUTO_TEST_CASE(ClampsOutsideTheRecordedRange)
{
    const auto data = makeRecordSeries();

    BOOST_CHECK_EQUAL(Opm::FluxBoundary::selectRecordAt(data, -5.0)->startTime, 0.0);
    BOOST_CHECK_EQUAL(Opm::FluxBoundary::selectRecordAt(data, 0.0)->startTime, 0.0);
    BOOST_CHECK_EQUAL(Opm::FluxBoundary::selectRecordAt(data, 1.0e6)->startTime, 30.0);
}

BOOST_AUTO_TEST_CASE(SelectsFromASingleRecord)
{
    auto data = makeRecordSeries();
    data.reportSteps.resize(1);

    BOOST_CHECK_EQUAL(Opm::FluxBoundary::selectRecordAt(data, -1.0)->startTime, 0.0);
    BOOST_CHECK_EQUAL(Opm::FluxBoundary::selectRecordAt(data, 5.0)->startTime, 0.0);
    BOOST_CHECK_EQUAL(Opm::FluxBoundary::selectRecordAt(data, 500.0)->startTime, 0.0);
}

BOOST_AUTO_TEST_CASE(SelectsNothingWithoutRecords)
{
    const Opm::EclIO::FluxFile::Data data;

    BOOST_CHECK(Opm::FluxBoundary::selectRecordAt(data, 1.0) == nullptr);
}

BOOST_AUTO_TEST_CASE(RecordsTileTheTimeAxisWithoutGaps)
{
    const auto data = makeRecordSeries();

    // Consecutive records must join exactly: a gap or an overlap would either
    // drop or double count the flow across the boundary.
    for (std::size_t i = 1; i < data.reportSteps.size(); ++i) {
        const auto previousEnd = data.reportSteps[i - 1].startTime
            + data.reportSteps[i - 1].stepLength;

        BOOST_CHECK_CLOSE(data.reportSteps[i].startTime, previousEnd, 1e-12);
    }
}
