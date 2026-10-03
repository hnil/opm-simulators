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

#include <config.h>

#define BOOST_TEST_MODULE ParentSummary

#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/flux/ParentSummary.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace {

// WOPR is a Rate, WOPT a Total, WBHP a Pressure.  Three samples at t = 10, 20
// and 30 seconds.  The cumulative is built from a constant 2.0 rate over the
// second interval so that linear interpolation has an exactly predictable
// midpoint.
Opm::EclIO::FluxFile::Data makeData()
{
    Opm::EclIO::FluxFile::Data data;

    data.summaryKeys = {"WOPR:P1", "WOPT:P1", "WBHP:P1"};
    data.summarySamples = {
        {10.0, {1.0, 100.0, 250.0}},
        {20.0, {2.0, 120.0, 260.0}},
        {30.0, {4.0, 160.0, 280.0}},
    };

    data.header.numSummaryKeys = 3;
    data.header.numSummarySamples = 3;
    data.header.summaryPerTimestep = true;

    return data;
}

Opm::ParentSummary makeSummary()
{
    auto summary = Opm::ParentSummary::fromFluxFile(makeData(), "TEST.FLUX");
    BOOST_REQUIRE(summary.has_value());
    return std::move(*summary);
}

} // namespace

BOOST_AUTO_TEST_CASE(LoadsKeysAndTimesFromFluxFile)
{
    const auto summary = makeSummary();

    BOOST_CHECK_EQUAL(summary.keys().size(), 3U);
    BOOST_CHECK_EQUAL(summary.times().size(), 3U);
    BOOST_CHECK_EQUAL(summary.sourceDescription(), "TEST.FLUX");
    BOOST_CHECK(!summary.empty());

    BOOST_CHECK(summary.has("WOPR:P1"));
    BOOST_CHECK(!summary.has("WGPR:P1"));
}

BOOST_AUTO_TEST_CASE(ReturnsNothingWithoutSummaryData)
{
    Opm::EclIO::FluxFile::Data data;
    BOOST_CHECK(!Opm::ParentSummary::fromFluxFile(data, "EMPTY.FLUX").has_value());

    auto keysOnly = makeData();
    keysOnly.summarySamples.clear();
    BOOST_CHECK(!Opm::ParentSummary::fromFluxFile(keysOnly, "EMPTY.FLUX").has_value());
}

BOOST_AUTO_TEST_CASE(ClassifiesKeywordTypes)
{
    using Type = Opm::SummaryConfigNode::Type;

    BOOST_CHECK(Opm::ParentSummary::isHeldConstant(Type::Rate));
    BOOST_CHECK(Opm::ParentSummary::isHeldConstant(Type::Mode));
    BOOST_CHECK(Opm::ParentSummary::isHeldConstant(Type::Count));

    BOOST_CHECK(!Opm::ParentSummary::isHeldConstant(Type::Total));
    BOOST_CHECK(!Opm::ParentSummary::isHeldConstant(Type::Pressure));

    // The keyword part of a composite key drives the classification.
    BOOST_CHECK(Opm::ParentSummary::keyType("WOPR:P1") == Type::Rate);
    BOOST_CHECK(Opm::ParentSummary::keyType("WOPT:P1") == Type::Total);
    BOOST_CHECK(Opm::ParentSummary::keyType("FOPR") == Type::Rate);
}

BOOST_AUTO_TEST_CASE(RateIsHeldPiecewiseConstant)
{
    const auto summary = makeSummary();

    // Exactly on a sample: that sample's value.
    BOOST_CHECK_CLOSE(summary.valueAt("WOPR:P1", 20.0), 2.0, 1e-12);
    BOOST_CHECK_CLOSE(summary.valueAt("WOPR:P1", 30.0), 4.0, 1e-12);

    // Strictly between samples: the LATER sample, because sample times are
    // interval end times and the rate stored at t[i] is the average over
    // (t[i-1], t[i]].  It must never be a blend of the two.
    BOOST_CHECK_CLOSE(summary.valueAt("WOPR:P1", 10.001), 2.0, 1e-12);
    BOOST_CHECK_CLOSE(summary.valueAt("WOPR:P1", 15.0), 2.0, 1e-12);
    BOOST_CHECK_CLOSE(summary.valueAt("WOPR:P1", 19.999), 2.0, 1e-12);
    BOOST_CHECK_CLOSE(summary.valueAt("WOPR:P1", 25.0), 4.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(RateHoldIsIntegralPreserving)
{
    const auto summary = makeSummary();

    // Holding the stored average across its own interval must reproduce the
    // parent's production over that interval.  Over (20, 30] the stored rate is
    // 4.0, and the cumulative grows by 160 - 120 = 40 over 10 seconds.
    const auto held = summary.valueAt("WOPR:P1", 25.0);
    const auto producedFromCumulative =
        summary.valueAt("WOPT:P1", 30.0) - summary.valueAt("WOPT:P1", 20.0);

    BOOST_CHECK_CLOSE(held * (30.0 - 20.0), producedFromCumulative, 1e-12);
}

BOOST_AUTO_TEST_CASE(RateOverAnIntervalIsItsAverage)
{
    const auto summary = makeSummary();

    // A step spanning several samples has to see what the parent produced
    // over the whole of it, not the rate of its last sample: that is what a
    // cumulative built from the rate adds up to. Over (10, 30] the parent made
    // 160 - 100 = 60 in 20 seconds, where the rate at 30 alone says 4.0.
    BOOST_CHECK_CLOSE(summary.valueOver("WOPR:P1", 10.0, 30.0), 3.0, 1e-12);
    BOOST_CHECK_CLOSE(summary.valueOver("WOPR:P1", 10.0, 30.0) * 20.0,
                      summary.valueAt("WOPT:P1", 30.0) - summary.valueAt("WOPT:P1", 10.0),
                      1e-12);

    // Straddling a sample, half in each interval.
    BOOST_CHECK_CLOSE(summary.valueOver("WOPR:P1", 15.0, 25.0), 3.0, 1e-12);

    // Inside a single interval it is that interval's value.
    BOOST_CHECK_CLOSE(summary.valueOver("WOPR:P1", 21.0, 29.0), 4.0, 1e-12);

    // Reaching past the last sample holds the last value there.
    BOOST_CHECK_CLOSE(summary.valueOver("WOPR:P1", 20.0, 40.0), 4.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(NonRateOverAnIntervalIsItsEndValue)
{
    const auto summary = makeSummary();

    BOOST_CHECK_CLOSE(summary.valueOver("WBHP:P1", 15.0, 25.0),
                      summary.valueAt("WBHP:P1", 25.0), 1e-12);
    BOOST_CHECK_CLOSE(summary.valueOver("WOPT:P1", 15.0, 25.0),
                      summary.valueAt("WOPT:P1", 25.0), 1e-12);

    // An empty interval is just a point.
    BOOST_CHECK_CLOSE(summary.valueOver("WOPR:P1", 25.0, 25.0),
                      summary.valueAt("WOPR:P1", 25.0), 1e-12);
}

BOOST_AUTO_TEST_CASE(TotalIsInterpolatedLinearly)
{
    const auto summary = makeSummary();

    BOOST_CHECK_CLOSE(summary.valueAt("WOPT:P1", 20.0), 120.0, 1e-12);

    // Midway between 120 and 160.
    BOOST_CHECK_CLOSE(summary.valueAt("WOPT:P1", 25.0), 140.0, 1e-12);

    // A quarter of the way between 100 and 120.
    BOOST_CHECK_CLOSE(summary.valueAt("WOPT:P1", 12.5), 105.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(PressureIsInterpolatedLinearly)
{
    const auto summary = makeSummary();

    BOOST_CHECK_CLOSE(summary.valueAt("WBHP:P1", 15.0), 255.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(ClampsOutsideTheSampledRange)
{
    const auto summary = makeSummary();

    // Before the first sample.
    BOOST_CHECK_CLOSE(summary.valueAt("WOPR:P1", 0.0), 1.0, 1e-12);
    BOOST_CHECK_CLOSE(summary.valueAt("WOPT:P1", -5.0), 100.0, 1e-12);

    // After the last sample.
    BOOST_CHECK_CLOSE(summary.valueAt("WOPR:P1", 1000.0), 4.0, 1e-12);
    BOOST_CHECK_CLOSE(summary.valueAt("WOPT:P1", 1000.0), 160.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(MissingKeyYieldsNaN)
{
    const auto summary = makeSummary();

    BOOST_CHECK(std::isnan(summary.valueAt("WGPR:P1", 15.0)));
}

BOOST_AUTO_TEST_CASE(MissingSummaryFileYieldsNothing)
{
    BOOST_CHECK(!Opm::ParentSummary::fromSummaryFile("no_such_file.SMSPEC").has_value());
}
