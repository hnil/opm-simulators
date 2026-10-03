/*
   Copyright 2026 Equinor ASA.

   This file is part of the Open Porous Media project (OPM).

   OPM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OPM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/io/eclipse/FluxFile.hpp>
#include <opm/io/eclipse/EclFile.hpp>

#include <algorithm>
#include <cstddef>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

namespace {

int fail(const std::string& msg)
{
    std::cerr << "flux_smoke_shape_check: " << msg << '\n';
    return EXIT_FAILURE;
}

bool parseInt(const std::string& text, int& value)
{
    std::istringstream is{text};
    is >> value;
    return is && is.eof();
}

std::string toUpper(std::string text)
{
    for (char& c : text) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return text;
}

bool parseMode(const std::string& text, Opm::EclIO::FluxFile::Mode& mode)
{
    const auto upper = toUpper(text);
    if (upper == "FLUX") {
        mode = Opm::EclIO::FluxFile::Mode::Flux;
        return true;
    }
    if (upper == "PRESSURE") {
        mode = Opm::EclIO::FluxFile::Mode::Pressure;
        return true;
    }
    if (upper == "BOTH") {
        mode = Opm::EclIO::FluxFile::Mode::Both;
        return true;
    }

    return false;
}

bool hasFluxMode(const Opm::EclIO::FluxFile::Mode mode)
{
    return (static_cast<int>(mode) & static_cast<int>(Opm::EclIO::FluxFile::Mode::Flux)) != 0;
}

bool hasPressureMode(const Opm::EclIO::FluxFile::Mode mode)
{
    return (static_cast<int>(mode) & static_cast<int>(Opm::EclIO::FluxFile::Mode::Pressure)) != 0;
}

} // namespace

int main(int argc, char** argv)
{
    // Arguments come in groups of 8 values:
    // <file.FLUX> <nx> <ny> <nz> <i1> <j1> <k1> <mode>
    constexpr int groupSize = 8;
    const int numArgs = argc - 1;

    if (numArgs < groupSize || (numArgs % groupSize) != 0) {
        return fail("usage: flux_smoke_shape_check <file.FLUX> <nx> <ny> <nz> <i1> <j1> <k1> <mode> [repeat]");
    }

    for (int base = 1; base < argc; base += groupSize) {
        const std::string path = argv[base + 0];

        int expectedNx = 0;
        int expectedNy = 0;
        int expectedNz = 0;
        int expectedI1 = 0;
        int expectedJ1 = 0;
        int expectedK1 = 0;
        Opm::EclIO::FluxFile::Mode expectedMode{};

        if (!parseInt(argv[base + 1], expectedNx)
            || !parseInt(argv[base + 2], expectedNy)
            || !parseInt(argv[base + 3], expectedNz)
            || !parseInt(argv[base + 4], expectedI1)
            || !parseInt(argv[base + 5], expectedJ1)
            || !parseInt(argv[base + 6], expectedK1)) {
            return fail("failed to parse integer arguments for " + path);
        }
        if (!parseMode(argv[base + 7], expectedMode)) {
            return fail("failed to parse mode argument for " + path);
        }

        const auto data = Opm::EclIO::FluxFile::read(path);
        const auto rawFile = Opm::EclIO::EclFile(path, Opm::EclIO::EclFile::Formatted{false}, true);
        const auto& h = data.header;

        if (h.parentNx != expectedNx || h.parentNy != expectedNy || h.parentNz != expectedNz) {
            return fail("unexpected parent dimensions in " + path);
        }

        if (h.boxI1 != expectedI1 || h.boxJ1 != expectedJ1 || h.boxK1 != expectedK1) {
            return fail("unexpected FLUXREG box start in " + path);
        }

        if (h.mode != expectedMode) {
            return fail("unexpected FLUXTYPE mode in " + path);
        }

        if (h.boxNx != 1 || h.boxNy != 1 || h.boxNz != 1) {
            return fail("unexpected FLUXREG box dimensions in " + path);
        }

        if (h.numCells != 1 || data.localToGlobal.size() != 1U) {
            return fail("unexpected cell mapping size in " + path);
        }

        if (h.numBoundaryFaces <= 0) {
            return fail("missing boundary faces in " + path);
        }

        if (data.boundaryFaces.size() != static_cast<std::size_t>(h.numBoundaryFaces)) {
            return fail("boundary face count mismatch in " + path);
        }

        if (h.numReportSteps <= 0) {
            return fail("missing report steps in " + path);
        }

        if (data.reportSteps.size() != static_cast<std::size_t>(h.numReportSteps)) {
            return fail("report step count mismatch in " + path);
        }

        const std::size_t expectedRates = static_cast<std::size_t>(h.numBoundaryFaces)
                                        * static_cast<std::size_t>(h.numPhases);

        const auto expectFlux = hasFluxMode(h.mode);
        const auto expectPressure = hasPressureMode(h.mode);

        // Surface-volume rates were dropped in format version 6. A flux-mode
        // file now carries component masses and nothing else, so that a
        // consumer cannot pick the one of the two that suits it and silently
        // disagree with the producer.
        if (rawFile.hasKey("FLXRATE")) {
            return fail("FLXRATE is no longer part of the format but is present in " + path);
        }
        if (rawFile.hasKey("FLXMASS") != expectFlux) {
            return fail("unexpected FLXMASS array presence in " + path);
        }
        if (rawFile.hasKey("FLXPRES") != expectPressure) {
            return fail("unexpected FLXPRES array presence in " + path);
        }
        if (rawFile.hasKey("FLXSATW") != (expectPressure && h.hasPhase(Opm::EclIO::FluxFile::Phase::Water))) {
            return fail("unexpected FLXSATW array presence in " + path);
        }
        if (rawFile.hasKey("FLXSATG") != (expectPressure && h.hasPhase(Opm::EclIO::FluxFile::Phase::Gas))) {
            return fail("unexpected FLXSATG array presence in " + path);
        }
        if (rawFile.hasKey("FLXRS") != expectPressure) {
            return fail("unexpected FLXRS array presence in " + path);
        }
        if (rawFile.hasKey("FLXRV") != expectPressure) {
            return fail("unexpected FLXRV array presence in " + path);
        }
        if (rawFile.hasKey("FLXTEMP") != h.hasTemperature) {
            return fail("unexpected FLXTEMP array presence in " + path);
        }

        for (const auto& step : data.reportSteps) {
            const auto expectedRateSize = expectFlux ? expectedRates : 0U;
            if (step.massRates.size() != expectedRateSize) {
                return fail("unexpected FLXMASS payload size in " + path);
            }

            const auto expectedFaceSize = expectPressure
                ? static_cast<std::size_t>(h.numBoundaryFaces)
                : 0U;

            if (step.pressures.size() != expectedFaceSize) {
                return fail("unexpected FLXPRES payload size in " + path);
            }
            if (step.swat.size() != expectedFaceSize) {
                return fail("unexpected FLXSATW payload size in " + path);
            }
            if (step.sgas.size() != expectedFaceSize) {
                return fail("unexpected FLXSATG payload size in " + path);
            }
            if (step.rs.size() != expectedFaceSize) {
                return fail("unexpected FLXRS payload size in " + path);
            }
            if (step.rv.size() != expectedFaceSize) {
                return fail("unexpected FLXRV payload size in " + path);
            }
        }

        if (expectFlux) {
            // A producer that fails to turn the FLORES computation on writes a
            // structurally valid file whose rates are all zero, which makes a
            // USEFLUX run behave as a closed region. Guard against that here
            // rather than leaving it to the equivalence tests.
            const auto anyRate =
                std::any_of(data.reportSteps.begin(), data.reportSteps.end(),
                            [](const auto& step)
                            {
                                return std::any_of(step.massRates.begin(), step.massRates.end(),
                                                   [](const double rate)
                                                   { return rate != 0.0; });
                            });

            if (!anyRate) {
                return fail("all boundary rates are zero in " + path);
            }
        }
    }

    return EXIT_SUCCESS;
}
