/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  OPM is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along
  with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

// Print the boundary data a FLUX file holds for one cell, as a time series.
//
// A FLUX file is stored face-major, which is the right layout for a consumer
// imposing boundary conditions and the wrong one for anybody trying to work
// out why a particular cell behaves the way it does. This turns it back round:
// name a cell and get every face on it, one row per record, in a shape that
// plots or diffs directly.

#include "config.h"

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/io/eclipse/FluxFile.hpp>

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;

constexpr double secondsPerDay = 24.0 * 60.0 * 60.0;

struct Options {
    std::string file;
    int i = 0;
    int j = 0;
    int k = 0;
    int precision = 6;
    std::string startDate;
    std::vector<std::string> summaryKeys;
    bool help = false;
};

void printUsage()
{
    std::cout
        << "usage: inspect_flux <file.FLUX> <i> <j> <k> [options]\n"
        << "       inspect_flux <file.FLUX> --summary=<KEY[,KEY...]> [options]\n"
        << "\n"
        << "Writes the boundary data a FLUX file carries for one cell as a time\n"
        << "series, one row per record in the file. The cell is given as a 1-based\n"
        << "cartesian triplet in the grid of the REDUCED run, which has the same\n"
        << "dimensions as the parent; only its active set is smaller.\n"
        << "\n"
        << "Columns are grouped by quantity with the faces cycling fastest, so a\n"
        << "cell with an I+ and a K- face gives I+:SWAT K-:SWAT I+:SGAS K-:SGAS and\n"
        << "so on. A face reached through a non-neighbour connection is labelled\n"
        << "NNC, numbered when the cell has more than one.\n"
        << "\n"
        << "Values are written in SI, as they are stored. The face geometry and\n"
        << "transmissibilities go to stderr, so that stdout stays a clean table.\n"
        << "\n"
        << "options:\n"
        << "  --summary=<KEY[,KEY...]>   Instead of a cell, list the embedded summary\n"
        << "                             samples of these keys, one row per sample,\n"
        << "                             in the parent's output units. The summary has\n"
        << "                             its own times, independent of the records.\n"
        << "  --start-date=<YYYY-MM-DD>  Report the first column as a date rather\n"
        << "                             than as days elapsed. A FLUX file records\n"
        << "                             no start date, so it has to be supplied.\n"
        << "  --precision=<n>            Significant digits (default 6)\n"
        << "  --help                     Show this message\n";
}

bool startsWith(const std::string& text, const std::string& prefix)
{
    return text.rfind(prefix, 0) == 0;
}

std::string valueAfterEquals(const std::string& arg)
{
    const auto pos = arg.find('=');
    if (pos == std::string::npos || pos + 1 >= arg.size()) {
        throw std::invalid_argument("expected --key=value form for argument '" + arg + "'");
    }
    return arg.substr(pos + 1);
}

int parsePositiveInt(const std::string& text, const char* what)
{
    try {
        std::size_t consumed = 0;
        const auto value = std::stoi(text, &consumed);
        if (consumed != text.size() || value <= 0) {
            throw std::invalid_argument("not a positive integer");
        }
        return value;
    }
    catch (const std::exception&) {
        throw std::invalid_argument(std::string{what} + " must be a positive integer, got '" + text + "'");
    }
}

Options parseOptions(int argc, char** argv)
{
    Options opt;
    std::vector<std::string> positional;

    for (int n = 1; n < argc; ++n) {
        const std::string arg = argv[n];

        if (arg == "--help" || arg == "-h") {
            opt.help = true;
        }
        else if (startsWith(arg, "--precision=")) {
            opt.precision = parsePositiveInt(valueAfterEquals(arg), "--precision");
        }
        else if (startsWith(arg, "--start-date=")) {
            opt.startDate = valueAfterEquals(arg);
        }
        else if (startsWith(arg, "--summary=")) {
            std::string keys = valueAfterEquals(arg);
            std::size_t pos = 0;
            while (pos <= keys.size()) {
                const auto comma = std::min(keys.find(',', pos), keys.size());
                if (comma > pos) {
                    opt.summaryKeys.push_back(keys.substr(pos, comma - pos));
                }
                pos = comma + 1;
            }
        }
        else if (startsWith(arg, "-")) {
            throw std::invalid_argument("unrecognized argument '" + arg + "'");
        }
        else {
            positional.push_back(arg);
        }
    }

    if (opt.help) {
        return opt;
    }

    if (!opt.summaryKeys.empty()) {
        if (positional.size() != 1U) {
            throw std::invalid_argument("--summary takes just <file.FLUX>, got "
                                        + std::to_string(positional.size()) + " positional argument(s)");
        }
        opt.file = positional[0];
        return opt;
    }

    if (positional.size() != 4U) {
        throw std::invalid_argument("expected <file.FLUX> <i> <j> <k>, got "
                                    + std::to_string(positional.size()) + " positional argument(s)");
    }

    opt.file = positional[0];
    opt.i = parsePositiveInt(positional[1], "i");
    opt.j = parsePositiveInt(positional[2], "j");
    opt.k = parsePositiveInt(positional[3], "k");

    return opt;
}

// Days since 1970-01-01, by Howard Hinnant's civil calendar algorithm. Used
// only to turn an elapsed time into a date when the caller supplies the start.
long long daysFromCivil(long long y, unsigned m, unsigned d)
{
    y -= m <= 2;
    const long long era = (y >= 0 ? y : y - 399) / 400;
    const auto yoe = static_cast<unsigned>(y - era * 400);
    const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;

    return era * 146097LL + static_cast<long long>(doe) - 719468LL;
}

std::string civilFromDays(long long z)
{
    z += 719468;
    const long long era = (z >= 0 ? z : z - 146096) / 146097;
    const auto doe = static_cast<unsigned long long>(z - era * 146097);
    const unsigned long long yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    const long long y = static_cast<long long>(yoe) + era * 400;
    const unsigned long long doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    const unsigned long long mp = (5 * doy + 2) / 153;
    const unsigned long long d = doy - (153 * mp + 2) / 5 + 1;
    const unsigned long long m = mp + (mp < 10 ? 3 : -9);

    std::ostringstream os;
    os << std::setfill('0') << std::setw(4) << (y + (m <= 2))
       << '-' << std::setw(2) << m
       << '-' << std::setw(2) << d;

    return os.str();
}

long long parseStartDate(const std::string& text)
{
    int year = 0;
    unsigned month = 0;
    unsigned day = 0;
    char dash1 = 0;
    char dash2 = 0;

    std::istringstream is(text);
    is >> year >> dash1 >> month >> dash2 >> day;
    if (!is || dash1 != '-' || dash2 != '-'
        || month < 1 || month > 12 || day < 1 || day > 31) {
        throw std::invalid_argument("--start-date must look like YYYY-MM-DD, got '" + text + "'");
    }

    return daysFromCivil(year, month, day);
}

std::string directionLabel(const int direction)
{
    switch (direction) {
    case static_cast<int>(Opm::FaceDir::XPlus):  return "I+";
    case static_cast<int>(Opm::FaceDir::XMinus): return "I-";
    case static_cast<int>(Opm::FaceDir::YPlus):  return "J+";
    case static_cast<int>(Opm::FaceDir::YMinus): return "J-";
    case static_cast<int>(Opm::FaceDir::ZPlus):  return "K+";
    case static_cast<int>(Opm::FaceDir::ZMinus): return "K-";
    default: return "NNC";
    }
}

std::string cartesianLabel(const int globalCell, const Opm::EclIO::FluxFile::Header& header)
{
    if (globalCell < 0) {
        return "none";
    }

    const auto nx = static_cast<long long>(header.parentNx);
    const auto ny = static_cast<long long>(header.parentNy);
    const auto cell = static_cast<long long>(globalCell);

    const auto i = cell % nx;
    const auto j = (cell / nx) % ny;
    const auto k = cell / (nx * ny);

    std::ostringstream os;
    os << '(' << (i + 1) << ',' << (j + 1) << ',' << (k + 1) << ')';

    return os.str();
}

//! \brief One output column: where to find its value in a record.
struct Column {
    std::string name;
    std::string unit;
    std::vector<double> Opm::EclIO::FluxFile::ReportStep::* member;
    std::size_t index;
    std::size_t width;
};

std::string formatDouble(const double value, const int precision)
{
    std::ostringstream os;
    os << std::setprecision(precision) << value;
    return os.str();
}

int run(const Options& opt)
{
    if (!fs::exists(opt.file)) {
        throw std::invalid_argument("FLUX file not found: '" + opt.file + "'");
    }

    const auto data = Opm::EclIO::FluxFile::read(opt.file);
    const auto& header = data.header;

    if (!opt.summaryKeys.empty()) {
        std::vector<std::size_t> columns;
        for (const auto& key : opt.summaryKeys) {
            const auto it = std::find(data.summaryKeys.begin(), data.summaryKeys.end(), key);
            if (it == data.summaryKeys.end()) {
                throw std::invalid_argument("summary key '" + key + "' is not embedded in '"
                                            + opt.file + "'");
            }
            columns.push_back(static_cast<std::size_t>(std::distance(data.summaryKeys.begin(), it)));
        }

        std::cout << "#" << std::setw(13) << "TIME";
        for (const auto& key : opt.summaryKeys) {
            std::cout << ' ' << std::setw(14) << key;
        }
        std::cout << "\n#" << std::setw(13) << "DAYS" << '\n';

        std::cout << std::setprecision(opt.precision);
        for (const auto& sample : data.summarySamples) {
            std::cout << std::setw(14) << sample.time / secondsPerDay;
            for (const auto column : columns) {
                std::cout << ' ' << std::setw(14) << sample.values[column];
            }
            std::cout << '\n';
        }

        return EXIT_SUCCESS;
    }

    const auto nx = static_cast<long long>(header.parentNx);
    const auto ny = static_cast<long long>(header.parentNy);
    const auto nz = static_cast<long long>(header.parentNz);

    if (opt.i > nx || opt.j > ny || opt.k > nz) {
        std::ostringstream os;
        os << "cell (" << opt.i << ',' << opt.j << ',' << opt.k
           << ") lies outside the parent grid, which is " << nx << 'x' << ny << 'x' << nz;
        throw std::invalid_argument(os.str());
    }

    const auto globalCell = static_cast<int>((opt.k - 1) * nx * ny + (opt.j - 1) * nx + (opt.i - 1));

    // The region description is box-local, so the cell has to be found in the
    // map rather than computed from the box origin: a box cell outside the
    // region carries -1 and must not be mistaken for a member of it.
    int localCell = -1;
    for (std::size_t local = 0; local < data.localToGlobal.size(); ++local) {
        if (data.localToGlobal[local] == globalCell) {
            localCell = static_cast<int>(local);
            break;
        }
    }

    if (localCell < 0) {
        std::ostringstream os;
        os << "cell (" << opt.i << ',' << opt.j << ',' << opt.k
           << ") is not part of the sector this file describes; its box starts at ("
           << header.boxI1 << ',' << header.boxJ1 << ',' << header.boxK1
           << ") and is " << header.boxNx << 'x' << header.boxNy << 'x' << header.boxNz;
        throw std::invalid_argument(os.str());
    }

    // Face order follows the file, because that order is what indexes the
    // payload arrays.
    std::vector<std::size_t> faceIndices;
    for (std::size_t f = 0; f < data.boundaryFaces.size(); ++f) {
        if (data.boundaryFaces[f].interiorLocalCell == localCell) {
            faceIndices.push_back(f);
        }
    }

    std::cerr << "inspect_flux: " << opt.file << '\n'
              << "  parent case      " << (data.names.empty() ? std::string{"unknown"} : data.names.front()) << '\n'
              << "  region           " << (data.names.size() > 1 ? data.names[1] : std::string{"unknown"}) << '\n'
              << "  parent grid      " << nx << 'x' << ny << 'x' << nz << '\n'
              << "  cell             (" << opt.i << ',' << opt.j << ',' << opt.k
              << ")  global " << globalCell << ", box-local " << localCell << '\n'
              << "  records          " << data.reportSteps.size() << '\n'
              << "  faces on cell    " << faceIndices.size()
              << " of " << data.boundaryFaces.size() << " in the region\n";

    if (faceIndices.empty()) {
        std::cerr << "\n  This cell has no boundary face. It is interior to the sector, or every\n"
                  << "  neighbour it has outside the sector is inactive or off the parent grid,\n"
                  << "  which the reduced run treats as a no-flow boundary.\n";
        return EXIT_SUCCESS;
    }

    // NNC faces need distinguishing from each other, since a cell may hold
    // several and they carry no direction to tell them apart.
    const auto nncCount = std::count_if(faceIndices.begin(), faceIndices.end(),
                                        [&data](const std::size_t f)
                                        {
                                            return data.boundaryFaces[f].direction
                                                == static_cast<int>(Opm::FaceDir::Unknown);
                                        });

    std::vector<std::string> faceLabels;
    faceLabels.reserve(faceIndices.size());
    int nncSeen = 0;
    for (const auto f : faceIndices) {
        const auto& face = data.boundaryFaces[f];
        auto label = directionLabel(face.direction);
        if (label == "NNC" && nncCount > 1) {
            label += std::to_string(++nncSeen);
        }
        faceLabels.push_back(label);
    }

    std::cerr << '\n';
    for (std::size_t n = 0; n < faceIndices.size(); ++n) {
        const auto& face = data.boundaryFaces[faceIndices[n]];
        std::cerr << "  " << std::left << std::setw(6) << faceLabels[n] << std::right
                  << " exterior cell " << std::setw(9) << face.exteriorGlobalCell
                  << ' ' << std::setw(14) << cartesianLabel(face.exteriorGlobalCell, header)
                  << "  transmissibility " << formatDouble(face.transmissibility, opt.precision)
                  << " m3.Pa.s"
                  << "  pvt region " << (face.exteriorPvtRegion + 1) << '\n';
    }
    std::cerr << '\n';

    std::vector<std::string> phaseTags;
    if (header.hasPhase(Opm::EclIO::FluxFile::Phase::Oil))   { phaseTags.emplace_back("O"); }
    if (header.hasPhase(Opm::EclIO::FluxFile::Phase::Water)) { phaseTags.emplace_back("W"); }
    if (header.hasPhase(Opm::EclIO::FluxFile::Phase::Gas))   { phaseTags.emplace_back("G"); }

    const auto numPhases = phaseTags.size();
    const auto valueWidth = static_cast<std::size_t>(opt.precision) + 8;

    // Present only what the file actually carries, so that a pressure-mode file
    // does not grow a wall of empty mass-rate columns.
    const auto carries = [&data](std::vector<double> Opm::EclIO::FluxFile::ReportStep::* member)
    {
        return std::any_of(data.reportSteps.begin(), data.reportSteps.end(),
                           [member](const auto& step) { return !(step.*member).empty(); });
    };

    std::vector<Column> columns;

    const auto addPerFace = [&](std::vector<double> Opm::EclIO::FluxFile::ReportStep::* member,
                                const std::string& quantity,
                                const std::string& unit)
    {
        if (!carries(member)) {
            return;
        }

        for (std::size_t n = 0; n < faceIndices.size(); ++n) {
            const auto name = faceLabels[n] + ':' + quantity;
            columns.push_back(Column{name, unit, member, faceIndices[n],
                                     std::max({name.size(), unit.size(), valueWidth})});
        }
    };

    const auto addPerFacePhase = [&](std::vector<double> Opm::EclIO::FluxFile::ReportStep::* member,
                                     const std::string& quantity,
                                     const std::string& unit)
    {
        if (!carries(member) || (numPhases == 0)) {
            return;
        }

        for (std::size_t p = 0; p < numPhases; ++p) {
            for (std::size_t n = 0; n < faceIndices.size(); ++n) {
                const auto name = faceLabels[n] + ':' + quantity + phaseTags[p];
                columns.push_back(Column{name, unit, member, faceIndices[n] * numPhases + p,
                                         std::max({name.size(), unit.size(), valueWidth})});
            }
        }
    };

    using Step = Opm::EclIO::FluxFile::ReportStep;
    addPerFace(&Step::pressures,   "PRES", "PA");
    addPerFace(&Step::swat,        "SWAT", "FRACTION");
    addPerFace(&Step::sgas,        "SGAS", "FRACTION");
    addPerFace(&Step::rs,          "RS",   "SM3/SM3");
    addPerFace(&Step::rv,          "RV",   "SM3/SM3");
    addPerFace(&Step::temperature, "TEMP", "K");
    addPerFacePhase(&Step::massRates,   "MASS", "KG/S");
    addPerFacePhase(&Step::relPerm,     "KR",   "FRACTION");
    addPerFacePhase(&Step::capPressure, "PC",   "PA");

    const bool asDate = !opt.startDate.empty();
    const auto startDay = asDate ? parseStartDate(opt.startDate) : 0LL;
    const std::string timeName = asDate ? "DATE" : "TIME";
    const std::string timeUnit = asDate ? "YYYY-MM-DD" : "DAYS";
    const auto timeWidth = std::max({timeName.size(), timeUnit.size(), valueWidth});

    // The '#' sits inside the first column so that the headings stay lined up
    // with the numbers underneath them.
    const auto headerCell = [](const std::string& text, const std::size_t width)
    {
        std::ostringstream os;
        os << std::setw(static_cast<int>(width) - 1) << std::right << text;
        return os.str();
    };

    std::ostringstream names;
    std::ostringstream units;
    names << '#' << headerCell(timeName, timeWidth);
    units << '#' << headerCell(timeUnit, timeWidth);
    for (const auto& column : columns) {
        names << ' ' << std::setw(static_cast<int>(column.width)) << std::right << column.name;
        units << ' ' << std::setw(static_cast<int>(column.width)) << std::right << column.unit;
    }

    std::cout << names.str() << '\n' << units.str() << '\n';

    for (const auto& step : data.reportSteps) {
        // A record covers [startTime, startTime + stepLength] and its state
        // arrays are the exterior values at the END of that interval, so that
        // is the instant the row belongs to.
        const auto endSeconds = step.startTime + step.stepLength;
        const auto endDays = endSeconds / secondsPerDay;

        std::ostringstream row;
        if (asDate) {
            row << std::setw(static_cast<int>(timeWidth)) << std::right
                << civilFromDays(startDay + static_cast<long long>(endDays));
        }
        else {
            row << std::setw(static_cast<int>(timeWidth)) << std::right
                << formatDouble(endDays, opt.precision);
        }

        for (const auto& column : columns) {
            const auto& values = step.*(column.member);
            const auto text = (column.index < values.size())
                ? formatDouble(values[column.index], opt.precision)
                : std::string{"-"};

            row << ' ' << std::setw(static_cast<int>(column.width)) << std::right << text;
        }

        std::cout << row.str() << '\n';
    }

    return EXIT_SUCCESS;
}

} // Anonymous namespace

int main(int argc, char** argv)
{
    try {
        const auto opt = parseOptions(argc, argv);
        if (opt.help || (argc < 2)) {
            printUsage();
            return EXIT_SUCCESS;
        }

        return run(opt);
    }
    catch (const std::exception& e) {
        std::cerr << "inspect_flux: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
