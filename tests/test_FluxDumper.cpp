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

#define BOOST_TEST_MODULE OPM_test_FluxDumper

#include <boost/test/unit_test.hpp>

#include <opm/io/eclipse/FluxFile.hpp>
#include <opm/simulators/flow/flux/FluxDumper.hpp>
#include <opm/simulators/flow/flux/FluxRegions.hpp>

#include <array>
#include <filesystem>
#include <stdexcept>
#include <vector>

namespace {

int globalIndex(const std::array<int, 3>& dims, const int i, const int j, const int k)
{
    return Opm::FluxRegions::cartesianIndex(dims, i - 1, j - 1, k - 1);
}

} // namespace

BOOST_AUTO_TEST_CASE(BuildsFluxFileMetadataFromRegion)
{
    const std::array<int, 3> dims{4, 1, 1};
    std::vector<int> regionValues(dims[0] * dims[1] * dims[2], 0);
    regionValues[globalIndex(dims, 2, 1, 1)] = 3;
    regionValues[globalIndex(dims, 3, 1, 1)] = 3;

    const auto regions = Opm::FluxRegions::extract(dims, regionValues);
    BOOST_REQUIRE_EQUAL(regions.size(), 1U);

    const std::vector<double> trans{11.0, 22.0};
    Opm::FluxDumper dumper("PARENT", 3, dims, regions.front(),
                           Opm::EclIO::FluxFile::Mode::Both,
                           Opm::EclIO::FluxFile::Sampling::Averaged,
                           static_cast<int>(Opm::EclIO::FluxFile::Phase::Oil)
                             | static_cast<int>(Opm::EclIO::FluxFile::Phase::Water),
                           trans);

    const auto& data = dumper.data();
    BOOST_CHECK_EQUAL(data.header.parentNx, 4);
    BOOST_CHECK_EQUAL(data.header.boxI1, 2);
    BOOST_CHECK_EQUAL(data.header.boxNx, 2);
    BOOST_CHECK_EQUAL(data.header.numBoundaryFaces, 2);
    BOOST_CHECK_EQUAL(data.header.numPhases, 2);
    BOOST_CHECK_EQUAL(data.localToGlobal.size(), 2U);
    BOOST_CHECK_EQUAL(data.boundaryFaces[0].transmissibility, 11.0);
    BOOST_CHECK_EQUAL(data.boundaryFaces[1].transmissibility, 22.0);
}

BOOST_AUTO_TEST_CASE(ValidatesReportStepVectorSizes)
{
    const std::array<int, 3> dims{2, 1, 1};
    std::vector<int> regionValues(dims[0] * dims[1] * dims[2], 0);
    regionValues[globalIndex(dims, 1, 1, 1)] = 1;

    const auto regions = Opm::FluxRegions::extract(dims, regionValues);
    BOOST_REQUIRE_EQUAL(regions.size(), 1U);

    Opm::FluxDumper dumper("PARENT", 1, dims, regions.front(),
                           Opm::EclIO::FluxFile::Mode::Both,
                           Opm::EclIO::FluxFile::Sampling::Instant,
                           static_cast<int>(Opm::EclIO::FluxFile::Phase::Oil)
                             | static_cast<int>(Opm::EclIO::FluxFile::Phase::Water)
                             | static_cast<int>(Opm::EclIO::FluxFile::Phase::Gas));

    Opm::FluxDumper::ReportStepData good;
    good.reportStep = 1;
    good.simStep = 3;
    good.startTime = 10.0;
    good.stepLength = 2.0;
    good.massRates = {1.0, 2.0, 3.0};
    good.pressures = {100.0};
    good.swat = {0.2};
    good.sgas = {0.1};
    good.rs = {20.0};
    good.rv = {0.5};
    good.temperature = {340.0};

    BOOST_CHECK_NO_THROW(dumper.appendReportStep(good));
    BOOST_CHECK_EQUAL(dumper.data().header.numReportSteps, 1);
    BOOST_CHECK(dumper.data().header.hasTemperature);

    Opm::FluxDumper::ReportStepData bad = good;
    bad.massRates = {1.0, 2.0};
    BOOST_CHECK_THROW(dumper.appendReportStep(bad), std::invalid_argument);
}

  BOOST_AUTO_TEST_CASE(AggregateRatesSupportsAveragedAndInstantSampling)
  {
    const std::vector<std::vector<double>> snapshots{
      {2.0, 10.0},
      {6.0, 4.0},
    };

    const auto instant = Opm::FluxDumper::aggregateRates(
      Opm::EclIO::FluxFile::Sampling::Instant,
      snapshots,
      {});
    BOOST_REQUIRE_EQUAL(instant.size(), 2U);
    BOOST_CHECK_CLOSE(instant[0], 6.0, 1e-12);
    BOOST_CHECK_CLOSE(instant[1], 4.0, 1e-12);

    const auto averaged = Opm::FluxDumper::aggregateRates(
      Opm::EclIO::FluxFile::Sampling::Averaged,
      snapshots,
      {1.0, 3.0});
    BOOST_REQUIRE_EQUAL(averaged.size(), 2U);
    BOOST_CHECK_CLOSE(averaged[0], 5.0, 1e-12);
    BOOST_CHECK_CLOSE(averaged[1], 5.5, 1e-12);
  }

  BOOST_AUTO_TEST_CASE(AggregateRatesValidatesInputs)
  {
    BOOST_CHECK_THROW(
      Opm::FluxDumper::aggregateRates(
        Opm::EclIO::FluxFile::Sampling::Averaged,
        {},
        {}),
      std::invalid_argument);

    BOOST_CHECK_THROW(
      Opm::FluxDumper::aggregateRates(
        Opm::EclIO::FluxFile::Sampling::Averaged,
        {{1.0}, {1.0, 2.0}},
        {1.0, 1.0}),
      std::invalid_argument);

    BOOST_CHECK_THROW(
      Opm::FluxDumper::aggregateRates(
        Opm::EclIO::FluxFile::Sampling::Averaged,
        {{1.0}, {2.0}},
        {2.0}),
      std::invalid_argument);

    BOOST_CHECK_THROW(
      Opm::FluxDumper::aggregateRates(
        Opm::EclIO::FluxFile::Sampling::Averaged,
        {{1.0}, {2.0}},
        {0.0, 0.0}),
      std::invalid_argument);
  }

  BOOST_AUTO_TEST_CASE(WritesFluxFileThroughFluxFileIo)
  {
    const std::array<int, 3> dims{2, 1, 1};
    std::vector<int> regionValues(dims[0] * dims[1] * dims[2], 0);
    regionValues[globalIndex(dims, 1, 1, 1)] = 4;

    const auto regions = Opm::FluxRegions::extract(dims, regionValues);
    BOOST_REQUIRE_EQUAL(regions.size(), 1U);

    Opm::FluxDumper dumper("PARENT", 4, dims, regions.front(),
                 Opm::EclIO::FluxFile::Mode::Flux,
                 Opm::EclIO::FluxFile::Sampling::Instant,
                 static_cast<int>(Opm::EclIO::FluxFile::Phase::Oil));

    Opm::FluxDumper::ReportStepData step;
    step.reportStep = 0;
    step.simStep = 0;
    step.startTime = 0.0;
    step.stepLength = 1.0;
    step.massRates = {12.5};
    dumper.appendReportStep(step);

    const auto outPath = std::filesystem::path{"test_fluxdumper_roundtrip.FLUX"};
    dumper.flush(outPath.string(), false);

    const auto loaded = Opm::EclIO::FluxFile::read(outPath.string());
    BOOST_CHECK_EQUAL(loaded.header.parentNx, 2);
    BOOST_CHECK_EQUAL(loaded.header.numBoundaryFaces, 1);
    BOOST_CHECK_EQUAL(loaded.header.numReportSteps, 1);
    BOOST_REQUIRE_EQUAL(loaded.reportSteps.size(), 1U);
    BOOST_REQUIRE_EQUAL(loaded.reportSteps[0].massRates.size(), 1U);
    BOOST_CHECK_CLOSE(loaded.reportSteps[0].massRates[0], 12.5, 1e-12);

    std::filesystem::remove(outPath);
  }

BOOST_AUTO_TEST_CASE(BuildsZeroFluxStepWithExpectedSizes)
{
  const std::array<int, 3> dims{4, 1, 1};
  std::vector<int> regionValues(dims[0] * dims[1] * dims[2], 0);
  regionValues[globalIndex(dims, 2, 1, 1)] = 2;
  regionValues[globalIndex(dims, 3, 1, 1)] = 2;

  const auto regions = Opm::FluxRegions::extract(dims, regionValues);
  BOOST_REQUIRE_EQUAL(regions.size(), 1U);

  Opm::FluxDumper dumper("PARENT", 2, dims, regions.front(),
               Opm::EclIO::FluxFile::Mode::Both,
               Opm::EclIO::FluxFile::Sampling::Averaged,
               static_cast<int>(Opm::EclIO::FluxFile::Phase::Oil)
               | static_cast<int>(Opm::EclIO::FluxFile::Phase::Water));

  const auto step = dumper.makeZeroFluxStep(4, 8, 12.0, 3.0);
  BOOST_CHECK_EQUAL(step.reportStep, 4);
  BOOST_CHECK_EQUAL(step.simStep, 8);
  BOOST_CHECK_CLOSE(step.startTime, 12.0, 1e-12);
  BOOST_CHECK_CLOSE(step.stepLength, 3.0, 1e-12);
  BOOST_CHECK_EQUAL(step.massRates.size(), 4U);
  BOOST_CHECK_EQUAL(step.pressures.size(), 2U);
  BOOST_CHECK_EQUAL(step.swat.size(), 2U);
  BOOST_CHECK_EQUAL(step.sgas.size(), 2U);
  BOOST_CHECK_EQUAL(step.rs.size(), 2U);
  BOOST_CHECK_EQUAL(step.rv.size(), 2U);
}

BOOST_AUTO_TEST_CASE(WritesAndReadsEmbeddedSummaryValues)
{
  const std::array<int, 3> dims{2, 1, 1};
  std::vector<int> regionValues(dims[0] * dims[1] * dims[2], 0);
  regionValues[globalIndex(dims, 1, 1, 1)] = 9;

  const auto regions = Opm::FluxRegions::extract(dims, regionValues);
  BOOST_REQUIRE_EQUAL(regions.size(), 1U);

  Opm::FluxDumper dumper("PARENT", 9, dims, regions.front(),
               Opm::EclIO::FluxFile::Mode::Pressure,
               Opm::EclIO::FluxFile::Sampling::Instant,
               static_cast<int>(Opm::EclIO::FluxFile::Phase::Oil));
  dumper.setSummaryKeys({"TIME", "FOPT"});

  Opm::FluxDumper::ReportStepData step;
  step.reportStep = 1;
  step.simStep = 1;
  step.startTime = 0.0;
  step.stepLength = 1.0;
  step.pressures = {123.0};
  step.rs = {0.0};
  step.rv = {0.0};
  dumper.appendReportStep(step);

  // Three samples against a single report step: the summary series is sampled
  // independently of the report step sequence.
  dumper.setSummaryMinSampleInterval(86400.0);
  dumper.appendSummarySample(0.25, {1.0, 2.0});
  dumper.appendSummarySample(0.50, {3.0, 4.0});
  dumper.appendSummarySample(1.00, {5.0, 6.0});

  const auto outPath = std::filesystem::path{"test_fluxdumper_summary_roundtrip.FLUX"};
  dumper.flush(outPath.string(), false);

  const auto loaded = Opm::EclIO::FluxFile::read(outPath.string());
  BOOST_REQUIRE_EQUAL(loaded.summaryKeys.size(), 2U);
  BOOST_CHECK_EQUAL(loaded.summaryKeys[0], "TIME");
  BOOST_CHECK_EQUAL(loaded.summaryKeys[1], "FOPT");
  BOOST_REQUIRE_EQUAL(loaded.reportSteps.size(), 1U);

  BOOST_REQUIRE_EQUAL(loaded.summarySamples.size(), 3U);
  BOOST_CHECK_CLOSE(loaded.summarySamples[0].time, 0.25, 1e-12);
  BOOST_CHECK_CLOSE(loaded.summarySamples[0].values[0], 1.0, 1e-12);
  BOOST_CHECK_CLOSE(loaded.summarySamples[0].values[1], 2.0, 1e-12);
  BOOST_CHECK_CLOSE(loaded.summarySamples[2].time, 1.00, 1e-12);
  BOOST_CHECK_CLOSE(loaded.summarySamples[2].values[0], 5.0, 1e-12);
  BOOST_CHECK_CLOSE(loaded.summarySamples[2].values[1], 6.0, 1e-12);

  BOOST_CHECK(loaded.header.summaryPerTimestep);
  BOOST_CHECK_CLOSE(loaded.header.summaryMinSampleInterval, 86400.0, 1e-12);

  std::filesystem::remove(outPath);
}

BOOST_AUTO_TEST_CASE(RejectsSummarySampleWithWrongWidth)
{
  const std::array<int, 3> dims{2, 1, 1};
  std::vector<int> regionValues(dims[0] * dims[1] * dims[2], 0);
  regionValues[globalIndex(dims, 1, 1, 1)] = 9;

  const auto regions = Opm::FluxRegions::extract(dims, regionValues);
  BOOST_REQUIRE_EQUAL(regions.size(), 1U);

  Opm::FluxDumper dumper("PARENT", 9, dims, regions.front(),
               Opm::EclIO::FluxFile::Mode::Pressure,
               Opm::EclIO::FluxFile::Sampling::Instant,
               static_cast<int>(Opm::EclIO::FluxFile::Phase::Oil));
  dumper.setSummaryKeys({"TIME", "FOPT"});

  BOOST_CHECK_THROW(dumper.appendSummarySample(0.5, {1.0}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(FlattensFaceMajorRatesInCanonicalPhaseOrder)
{
  const std::array<int, 3> dims{4, 1, 1};
  std::vector<int> regionValues(dims[0] * dims[1] * dims[2], 0);
  regionValues[globalIndex(dims, 2, 1, 1)] = 6;
  regionValues[globalIndex(dims, 3, 1, 1)] = 6;

  const auto regions = Opm::FluxRegions::extract(dims, regionValues);
  BOOST_REQUIRE_EQUAL(regions.size(), 1U);

  Opm::FluxDumper dumper("PARENT", 6, dims, regions.front(),
               Opm::EclIO::FluxFile::Mode::Flux,
               Opm::EclIO::FluxFile::Sampling::Instant,
               static_cast<int>(Opm::EclIO::FluxFile::Phase::Oil)
               | static_cast<int>(Opm::EclIO::FluxFile::Phase::Gas));

  const auto rates = dumper.makeFaceMajorRates(
    [](const Opm::FluxRegions::BoundaryFace& face, const Opm::EclIO::FluxFile::Phase phase)
    {
      const int phaseCode = (phase == Opm::EclIO::FluxFile::Phase::Oil) ? 10 : 30;
      return static_cast<double>(face.interiorGlobalCell * 100 + phaseCode);
    });

  BOOST_REQUIRE_EQUAL(rates.size(), 4U);
  BOOST_CHECK_CLOSE(rates[0], globalIndex(dims, 2, 1, 1) * 100.0 + 10.0, 1e-12);
  BOOST_CHECK_CLOSE(rates[1], globalIndex(dims, 2, 1, 1) * 100.0 + 30.0, 1e-12);
  BOOST_CHECK_CLOSE(rates[2], globalIndex(dims, 3, 1, 1) * 100.0 + 10.0, 1e-12);
  BOOST_CHECK_CLOSE(rates[3], globalIndex(dims, 3, 1, 1) * 100.0 + 30.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(FlattensNncBoundaryFaceRates)
{
  const std::array<int, 3> dims{3, 1, 1};
  std::vector<int> regionValues(dims[0] * dims[1] * dims[2], 0);
  regionValues[globalIndex(dims, 1, 1, 1)] = 8;

  const std::vector<std::array<int, 2>> nnc{
    {globalIndex(dims, 1, 1, 1), globalIndex(dims, 3, 1, 1)}
  };

  const auto regions = Opm::FluxRegions::extract(dims, regionValues, nnc);
  BOOST_REQUIRE_EQUAL(regions.size(), 1U);

  Opm::FluxDumper dumper("PARENT", 8, dims, regions.front(),
               Opm::EclIO::FluxFile::Mode::Flux,
               Opm::EclIO::FluxFile::Sampling::Instant,
               static_cast<int>(Opm::EclIO::FluxFile::Phase::Oil));

  const auto rates = dumper.makeFaceMajorRates(
    [](const Opm::FluxRegions::BoundaryFace& face, const Opm::EclIO::FluxFile::Phase)
    {
      return face.isNnc ? 99.0 : 11.0;
    });

  BOOST_REQUIRE_EQUAL(rates.size(), 2U);
  BOOST_CHECK_CLOSE(rates[0], 11.0, 1e-12);
  BOOST_CHECK_CLOSE(rates[1], 99.0, 1e-12);
}
