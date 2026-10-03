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

#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace {

int fail(const std::string& msg)
{
    std::cerr << "flux_smoke_mutate: " << msg << '\n';
    return EXIT_FAILURE;
}

bool hasFluxMode(const Opm::EclIO::FluxFile::Mode mode)
{
    return (static_cast<int>(mode) & static_cast<int>(Opm::EclIO::FluxFile::Mode::Flux)) != 0;
}

// A FLUX file is only ever produced by appending, so even a caller that holds
// the whole thing in hand goes through the writer.
void writeAll(const std::string& path, const Opm::EclIO::FluxFile::Data& data)
{
    Opm::EclIO::FluxFile::Writer writer(path, /*formatted=*/false, data);
    writer.appendRecords(data.reportSteps);
    writer.appendSummarySamples(data.summarySamples);
    writer.close();
}

double mutationValue(const std::string& mutation)
{
    if (mutation == "nan") {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (mutation == "inf") {
        return std::numeric_limits<double>::infinity();
    }
    if (mutation == "ninf") {
        return -std::numeric_limits<double>::infinity();
    }

    throw std::invalid_argument("unknown mutation '" + mutation + "' (expected nan|inf|ninf)");
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 4) {
        return fail("usage: flux_smoke_mutate <input.FLUX> <output.FLUX> "
                    "<nan|inf|ninf|strip-summary>");
    }

    const std::string inputPath = argv[1];
    const std::string outputPath = argv[2];
    const std::string mutation = argv[3];

    auto data = Opm::EclIO::FluxFile::read(inputPath);

    if (mutation == "strip-summary") {
        // Remove the embedded parent summary vectors so that a USEFLUX run has
        // to fall back to the parent run's own summary output.
        if (data.summaryKeys.empty() && data.summarySamples.empty()) {
            return fail("input file already has no embedded summary data");
        }

        data.summaryKeys.clear();
        data.summarySamples.clear();
        data.header.numSummaryKeys = 0;
        data.header.numSummarySamples = 0;
        data.header.summaryPerTimestep = false;
        data.header.summaryMinSampleInterval = 0.0;

        writeAll(outputPath, data);
        return EXIT_SUCCESS;
    }

    if (!hasFluxMode(data.header.mode)) {
        return fail("input file does not contain FLUX mode rates");
    }

    if (data.reportSteps.empty()) {
        return fail("input file has no report steps");
    }

    if (data.reportSteps.front().massRates.empty()) {
        return fail("input file has empty FLXRATE payload");
    }

    try {
        data.reportSteps.front().massRates.front() = mutationValue(mutation);
    }
    catch (const std::exception& e) {
        return fail(e.what());
    }

    writeAll(outputPath, data);
    return EXIT_SUCCESS;
}
