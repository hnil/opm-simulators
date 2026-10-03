/*
   Copyright 2026 Equinor ASA.

   This file is part of the Open Porous Media project (OPM).

   OPM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
*/

#include <opm/io/eclipse/FluxFile.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

int fail(const std::string& msg)
{
    std::cerr << "flux_smoke_compare: " << msg << '\n';
    return EXIT_FAILURE;
}

bool closeEnough(const double lhs, const double rhs, const double tol)
{
    const auto scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) <= tol * scale;
}

//! \brief Compare two values on their own scale rather than against unity.
//!
//! \details Needed for quantities that are genuinely tiny in SI.
//!   Transmissibility is the one that matters here: it runs around 1e-12, so
//!   closeEnough() treats every value as indistinguishable from zero and from
//!   every other value. Returns true when both are zero, which is how a file
//!   that carries no transmissibility at all is represented.
bool relativelyClose(const double lhs, const double rhs, const double tol)
{
    const auto scale = std::max(std::abs(lhs), std::abs(rhs));
    return std::abs(lhs - rhs) <= tol * scale;
}

bool compareVector(const std::vector<double>& lhs,
                   const std::vector<double>& rhs,
                   const std::string& name,
                   const double tol)
{
    if (lhs.size() != rhs.size()) {
        std::cerr << name << " size mismatch: " << lhs.size() << " vs " << rhs.size() << '\n';
        return false;
    }

    for (std::size_t i = 0; i < lhs.size(); ++i) {
        if (!closeEnough(lhs[i], rhs[i], tol)) {
            std::cerr << name << " mismatch at index " << i << ": " << lhs[i] << " vs " << rhs[i] << '\n';
            return false;
        }
    }

    return true;
}

bool compareStep(const Opm::EclIO::FluxFile::ReportStep& lhs,
                 const Opm::EclIO::FluxFile::ReportStep& rhs,
                 const std::size_t stepIdx,
                 const bool ignorePressures,
                 const bool ignoreTimes)
{
    if (lhs.reportStep != rhs.reportStep || lhs.simStep != rhs.simStep) {
        std::cerr << "step meta mismatch at step " << stepIdx
                  << ": reportStep " << lhs.reportStep << " vs " << rhs.reportStep
                  << ", simStep " << lhs.simStep << " vs " << rhs.simStep << '\n';
        return false;
    }
    if (!ignoreTimes
        && (!closeEnough(lhs.startTime, rhs.startTime, 1e-10)
            || !closeEnough(lhs.stepLength, rhs.stepLength, 1e-10))) {
        std::cerr << "step time mismatch at step " << stepIdx
                  << ": startTime " << lhs.startTime << " vs " << rhs.startTime
                  << ", stepLength " << lhs.stepLength << " vs " << rhs.stepLength << '\n';
        return false;
    }

    const bool pressureMatch = ignorePressures
        ? lhs.pressures.size() == rhs.pressures.size()
        : compareVector(lhs.pressures, rhs.pressures, "pressures", 1e-2);
    if (ignorePressures && !pressureMatch) {
        std::cerr << "pressures size mismatch: " << lhs.pressures.size() << " vs " << rhs.pressures.size() << '\n';
        return false;
    }

    // Saturations get single precision's worth of room: a restart file holds
    // them as floats, so a file make_flux built from one agrees with DUMPFLUX's
    // doubles only to about 1e-7.
    const bool basicMatch = compareVector(lhs.massRates, rhs.massRates, "massRates", 1e-10)
        && pressureMatch
        && compareVector(lhs.swat, rhs.swat, "swat", 1e-6)
        && compareVector(lhs.sgas, rhs.sgas, "sgas", 1e-6)
        && compareVector(lhs.rs, rhs.rs, "rs", 1e-4)
        && compareVector(lhs.rv, rhs.rv, "rv", 1e-4)
        && compareVector(lhs.temperature, rhs.temperature, "temperature", 1e-6)
        && compareVector(lhs.relPerm, rhs.relPerm, "relPerm", 1e-5)
        && compareVector(lhs.capPressure, rhs.capPressure, "capPressure", 1e-2);

    return basicMatch;
}

} // namespace

int main(int argc, char** argv)
{
    bool ignorePressures = false;
    bool ignoreTimes = false;
    for (int i = 3; i < argc; ++i) {
        const std::string arg{argv[i]};
        if (arg == "--ignore-pressures") {
            ignorePressures = true;
            continue;
        }
        if (arg == "--ignore-times") {
            ignoreTimes = true;
            continue;
        }

        return fail("usage: flux_smoke_compare <expected.FLUX> <actual.FLUX>"
                    " [--ignore-pressures] [--ignore-times]");
    }

    if (argc < 3) {
        return fail("usage: flux_smoke_compare <expected.FLUX> <actual.FLUX>"
                    " [--ignore-pressures] [--ignore-times]");
    }

    const auto expected = Opm::EclIO::FluxFile::read(argv[1]);
    const auto actual = Opm::EclIO::FluxFile::read(argv[2]);

    std::size_t expectedStepOffset = 0;
    if (expected.reportSteps.size() == actual.reportSteps.size() + 1
        && !expected.reportSteps.empty()
        && expected.reportSteps.front().reportStep == 0) {
        expectedStepOffset = 1;
    }

    const bool headerMatch = expected.header.version == actual.header.version
        && expected.header.parentNx == actual.header.parentNx
        && expected.header.parentNy == actual.header.parentNy
        && expected.header.parentNz == actual.header.parentNz
        && expected.header.boxI1 == actual.header.boxI1
        && expected.header.boxJ1 == actual.header.boxJ1
        && expected.header.boxK1 == actual.header.boxK1
        && expected.header.boxNx == actual.header.boxNx
        && expected.header.boxNy == actual.header.boxNy
        && expected.header.boxNz == actual.header.boxNz
        && expected.header.numCells == actual.header.numCells
        && expected.header.numBoundaryFaces == actual.header.numBoundaryFaces
        && expected.header.numReportSteps == static_cast<int>(actual.reportSteps.size() + expectedStepOffset)
        && actual.header.numReportSteps == static_cast<int>(actual.reportSteps.size())
        && expected.header.numPhases == actual.header.numPhases
        && expected.header.phaseMask == actual.header.phaseMask
        && expected.header.hasTemperature == actual.header.hasTemperature
        && expected.header.mode == actual.header.mode
        && expected.header.sampling == actual.header.sampling;

    if (!headerMatch) {
        const auto& lhs = expected.header;
        const auto& rhs = actual.header;
        std::cerr << "header mismatch:\n"
                  << "  version: " << lhs.version << " vs " << rhs.version << '\n'
                  << "  parentNx: " << lhs.parentNx << " vs " << rhs.parentNx << '\n'
                  << "  parentNy: " << lhs.parentNy << " vs " << rhs.parentNy << '\n'
                  << "  parentNz: " << lhs.parentNz << " vs " << rhs.parentNz << '\n'
                  << "  boxI1: " << lhs.boxI1 << " vs " << rhs.boxI1 << '\n'
                  << "  boxJ1: " << lhs.boxJ1 << " vs " << rhs.boxJ1 << '\n'
                  << "  boxK1: " << lhs.boxK1 << " vs " << rhs.boxK1 << '\n'
                  << "  boxNx: " << lhs.boxNx << " vs " << rhs.boxNx << '\n'
                  << "  boxNy: " << lhs.boxNy << " vs " << rhs.boxNy << '\n'
                  << "  boxNz: " << lhs.boxNz << " vs " << rhs.boxNz << '\n'
                  << "  numCells: " << lhs.numCells << " vs " << rhs.numCells << '\n'
                  << "  numBoundaryFaces: " << lhs.numBoundaryFaces << " vs " << rhs.numBoundaryFaces << '\n'
                  << "  numReportSteps: " << lhs.numReportSteps << " vs " << rhs.numReportSteps << '\n'
                  << "  numPhases: " << lhs.numPhases << " vs " << rhs.numPhases << '\n'
                  << "  phaseMask: " << lhs.phaseMask << " vs " << rhs.phaseMask << '\n'
                  << "  hasTemperature: " << lhs.hasTemperature << " vs " << rhs.hasTemperature << '\n'
                  << "  mode: " << static_cast<int>(lhs.mode) << " vs " << static_cast<int>(rhs.mode) << '\n'
                  << "  sampling: " << static_cast<int>(lhs.sampling) << " vs " << static_cast<int>(rhs.sampling) << '\n';
        return fail("header mismatch");
    }
    if (expected.names != actual.names) {
        return fail("name metadata mismatch");
    }
    if (expected.localToGlobal != actual.localToGlobal) {
        return fail("LOC2GLOB mismatch");
    }
    const bool compareSummary = !expected.summaryKeys.empty();
    if (compareSummary && expected.summaryKeys != actual.summaryKeys) {
        return fail("summary key mismatch");
    }
    if (expected.boundaryFaces.size() != actual.boundaryFaces.size()) {
        return fail("boundary face count mismatch");
    }

    for (std::size_t i = 0; i < expected.boundaryFaces.size(); ++i) {
        const auto& lhs = expected.boundaryFaces[i];
        const auto& rhs = actual.boundaryFaces[i];
        // A tool that rebuilds the boundary from a parent's INIT file reads
        // single-precision transmissibilities where the simulator had doubles,
        // so allow for that but nothing looser.
        if (lhs.interiorLocalCell != rhs.interiorLocalCell
            || lhs.direction != rhs.direction
            || lhs.exteriorGlobalCell != rhs.exteriorGlobalCell
            || !relativelyClose(lhs.transmissibility, rhs.transmissibility, 1e-6)) {
            std::cerr << "boundary face mismatch at index " << i << ":\n"
                      << "  interiorLocalCell: " << lhs.interiorLocalCell << " vs " << rhs.interiorLocalCell << '\n'
                      << "  direction: " << lhs.direction << " vs " << rhs.direction << '\n'
                      << "  exteriorGlobalCell: " << lhs.exteriorGlobalCell << " vs " << rhs.exteriorGlobalCell << '\n'
                      << "  transmissibility: " << lhs.transmissibility << " vs " << rhs.transmissibility << '\n';
            return fail("boundary face mismatch at index " + std::to_string(i));
        }
    }

    if (expected.reportSteps.size() != actual.reportSteps.size() + expectedStepOffset) {
        return fail("report step count mismatch");
    }

    for (std::size_t i = 0; i < expected.reportSteps.size(); ++i) {
        if (i < expectedStepOffset) {
            continue;
        }

        if (!compareStep(expected.reportSteps[i],
                         actual.reportSteps[i - expectedStepOffset],
                         i - expectedStepOffset,
                         ignorePressures,
                         ignoreTimes)) {
            return fail("report step payload mismatch at index " + std::to_string(i));
        }
    }

    // Summary samples form their own series, independent of the report steps.
    if (compareSummary) {
        if (expected.summarySamples.size() != actual.summarySamples.size()) {
            std::cerr << "summary sample count mismatch: "
                      << expected.summarySamples.size() << " vs "
                      << actual.summarySamples.size() << '\n';
            return fail("summary sample count mismatch");
        }

        for (std::size_t i = 0; i < expected.summarySamples.size(); ++i) {
            const auto& lhs = expected.summarySamples[i];
            const auto& rhs = actual.summarySamples[i];

            if (!ignoreTimes && !closeEnough(lhs.time, rhs.time, 1e-10)) {
                std::cerr << "summary sample time mismatch at index " << i << ": "
                          << lhs.time << " vs " << rhs.time << '\n';
                return fail("summary sample time mismatch");
            }

            // A summary file holds floats, as a restart file does.
            if (!compareVector(lhs.values, rhs.values, "summaryValues", 1e-6)) {
                return fail("summary sample payload mismatch at index " + std::to_string(i));
            }
        }
    }

    return EXIT_SUCCESS;
}
