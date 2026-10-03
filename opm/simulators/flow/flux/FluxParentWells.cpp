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

#include <opm/simulators/flow/flux/FluxParentWells.hpp>

#include <opm/common/utility/TimeService.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <fmt/format.h>

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <unordered_map>

namespace {

int wellTypeCode(const Opm::Well& well)
{
    if (well.isProducer()) {
        return 1;
    }

    switch (well.injectorType()) {
    case Opm::InjectorType::GAS: return 3;
    case Opm::InjectorType::OIL: return 4;
    default:                     return 2;
    }
}

int wellStatusCode(const Opm::Well::Status status)
{
    switch (status) {
    case Opm::Well::Status::OPEN: return 1;
    case Opm::Well::Status::STOP: return 2;
    case Opm::Well::Status::SHUT: return 3;
    default:                      return 4;
    }
}

const char* wellStatusName(const int code)
{
    switch (code) {
    case 1:  return "OPEN";
    case 2:  return "STOP";
    case 3:  return "SHUT";
    default: return "AUTO";
    }
}

const char* injectedPhaseName(const int type)
{
    switch (type) {
    case 3:  return "GAS";
    case 4:  return "OIL";
    default: return "WATER";
    }
}

std::vector<std::unique_ptr<Opm::DeckKeyword>> parseKeywords(const std::string& text)
{
    // Depths are written in metres whatever the deck's own units, so the
    // snippet says METRIC and the parser converts to SI from that.
    Opm::Parser parser;
    parser.silent(true);
    auto deck = parser.parseString("METRIC\n\n" + text);
    deck.remove_keywords(0, 1);

    std::vector<std::unique_ptr<Opm::DeckKeyword>> keywords;
    for (const auto& keyword : deck) {
        keywords.push_back(std::make_unique<Opm::DeckKeyword>(keyword));
    }

    return keywords;
}

// The producer's report step that is in force at time t: the last one that
// starts at or before it.
std::ptrdiff_t parentStepAt(const std::vector<double>& stepStart, const double t)
{
    std::ptrdiff_t step = -1;
    for (std::size_t s = 0; s < stepStart.size(); ++s) {
        if (stepStart[s] <= t + 0.5) {
            step = static_cast<std::ptrdiff_t>(s);
        }
    }

    return step;
}

} // namespace

namespace Opm {

EclIO::FluxFile::ParentWells describeFluxParentWells(const Schedule& schedule)
{
    EclIO::FluxFile::ParentWells parent;

    parent.wells = schedule.wellNames();
    parent.groups = schedule.groupNames();
    if (parent.wells.empty()) {
        return {};
    }

    std::unordered_map<std::string, int> groupIndex;
    for (std::size_t g = 0; g < parent.groups.size(); ++g) {
        groupIndex.emplace(parent.groups[g], static_cast<int>(g) + 1);
    }

    const auto indexOf = [&groupIndex](const std::string& name)
    {
        const auto it = groupIndex.find(name);
        return (it == groupIndex.end()) ? 0 : it->second;
    };

    const auto numSteps = schedule.size();
    parent.wellRefDepth.assign(parent.wells.size(), std::numeric_limits<double>::quiet_NaN());

    for (std::size_t step = 0; step < numSteps; ++step) {
        parent.stepStart.push_back(static_cast<double>(TimeService::to_time_t(schedule.simTime(step))));

        for (std::size_t w = 0; w < parent.wells.size(); ++w) {
            const auto& name = parent.wells[w];
            if (!schedule.hasWell(name, step)) {
                parent.wellGroup.push_back(0);
                parent.wellType.push_back(0);
                parent.wellStatus.push_back(0);
                parent.wellEfficiency.push_back(1.0);
                continue;
            }

            const auto& well = schedule.getWell(name, step);
            parent.wellGroup.push_back(indexOf(well.groupName()));
            parent.wellType.push_back(wellTypeCode(well));
            parent.wellStatus.push_back(wellStatusCode(well.getStatus()));
            parent.wellEfficiency.push_back(well.getEfficiencyFactor());

            if (well.hasRefDepth()) {
                parent.wellRefDepth[w] = well.getRefDepth();
            }
        }

        for (const auto& name : parent.groups) {
            if (!schedule.hasGroup(name, step)) {
                parent.groupParent.push_back(0);
                parent.groupEfficiency.push_back(1.0);
                continue;
            }

            const auto& group = schedule.getGroup(name, step);
            parent.groupParent.push_back((name == "FIELD") ? 0 : indexOf(group.parent()));
            parent.groupEfficiency.push_back(group.getGroupEfficiencyFactor());
        }
    }

    return parent;
}

std::vector<std::string>
addFluxParentWells(Schedule& schedule,
                   const EclIO::FluxFile::ParentWells& parent,
                   const ParseContext& parseContext,
                   ErrorGuard& errors)
{
    std::vector<std::string> added;
    if (parent.empty()) {
        return added;
    }

    const auto numWells = parent.wells.size();
    const auto numGroups = parent.groups.size();

    std::vector<std::size_t> missingWells;
    for (std::size_t w = 0; w < numWells; ++w) {
        if (!schedule.hasWell(parent.wells[w])) {
            missingWells.push_back(w);
            added.push_back(parent.wells[w]);
        }
    }

    if (missingWells.empty()) {
        return added;
    }

    const auto& phases = schedule.runspec().phases();
    const auto producerPhase = phases.active(Phase::OIL)
        ? "OIL" : (phases.active(Phase::GAS) ? "GAS" : "WATER");

    // What has been put into the schedule so far, so that a keyword goes in
    // only where something changes.
    struct WellState { int group = 0; int type = 0; int status = 0; double efficiency = 1.0; };
    std::vector<WellState> wellState(numWells);
    std::vector<int> groupParent(numGroups, -1);
    std::vector<double> groupEfficiency(numGroups, 1.0);

    std::map<std::size_t, std::vector<DeckKeyword>> keywordsByStep;

    for (std::size_t step = 0; step < schedule.size(); ++step) {
        const auto p = parentStepAt(parent.stepStart,
                                    static_cast<double>(TimeService::to_time_t(schedule.simTime(step))));
        if (p < 0) {
            continue;
        }

        const auto groupAt = [&parent, p, numGroups](const std::size_t g)
        { return parent.groupParent[static_cast<std::size_t>(p) * numGroups + g]; };

        std::string gruptree;
        std::string gefac;
        for (std::size_t g = 0; g < numGroups; ++g) {
            const auto& name = parent.groups[g];
            const auto parentIndex = groupAt(g);

            // FIELD is always there, and a group the reduced run defines
            // itself is left as it defined it.
            if ((name == "FIELD") || (parentIndex == 0)
                || ((groupParent[g] < 0) && schedule.hasGroup(name, step)))
            {
                continue;
            }

            if (groupParent[g] != parentIndex) {
                gruptree += fmt::format(" '{}' '{}' /\n", name,
                                        parent.groups[static_cast<std::size_t>(parentIndex) - 1]);
                groupParent[g] = parentIndex;
            }

            const auto efficiency = parent.groupEfficiency[static_cast<std::size_t>(p) * numGroups + g];
            if (efficiency != groupEfficiency[g]) {
                gefac += fmt::format(" '{}' {} /\n", name, efficiency);
                groupEfficiency[g] = efficiency;
            }
        }

        std::string welspecs;
        std::string wconprod;
        std::string wconinje;
        std::string wefac;
        for (const auto w : missingWells) {
            const auto slot = static_cast<std::size_t>(p) * numWells + w;
            const auto group = parent.wellGroup[slot];
            if (group == 0) {
                continue;
            }

            const auto& name = parent.wells[w];
            const auto type = parent.wellType[slot];
            const auto status = parent.wellStatus[slot];
            auto& state = wellState[w];

            if (state.group != group) {
                // The well head is put at (1,1): the well has no connections
                // here, and the parent's head need not lie in this grid.
                const auto depth = parent.wellRefDepth[w];
                welspecs += fmt::format(" '{}' '{}' 1 1 {} '{}' /\n", name,
                                        parent.groups[static_cast<std::size_t>(group) - 1],
                                        std::isfinite(depth) ? fmt::format("{}", depth) : std::string{"1*"},
                                        (type <= 1) ? producerPhase : injectedPhaseName(type));
            }

            if ((state.type != type) || (state.status != status)) {
                // The control is a placeholder. A well without connections is
                // never solved for; its rates come from the parent summary.
                if (type <= 1) {
                    wconprod += fmt::format(" '{}' '{}' 'ORAT' 0.0 /\n", name, wellStatusName(status));
                }
                else {
                    wconinje += fmt::format(" '{}' '{}' '{}' 'RATE' 0.0 /\n", name,
                                            injectedPhaseName(type), wellStatusName(status));
                }
            }

            const auto efficiency = parent.wellEfficiency[slot];
            if (efficiency != state.efficiency) {
                wefac += fmt::format(" '{}' {} /\n", name, efficiency);
            }

            state = WellState{group, type, status, efficiency};
        }

        std::string text;
        const auto append = [&text](const char* keyword, const std::string& records)
        {
            if (!records.empty()) {
                text += fmt::format("{}\n{}/\n\n", keyword, records);
            }
        };

        append("GRUPTREE", gruptree);
        append("GEFAC", gefac);
        append("WELSPECS", welspecs);
        append("WCONPROD", wconprod);
        append("WCONINJE", wconinje);
        append("WEFAC", wefac);

        if (!text.empty()) {
            auto& keywords = keywordsByStep[step];
            for (auto& keyword : parseKeywords(text)) {
                keywords.push_back(std::move(*keyword));
            }
        }
    }

    // All at once: the schedule is rebuilt a single time, and every keyword
    // of the deck is then handled with every added well in place. Inserted
    // step by step, each rebuild would meet the deck's references to wells not
    // yet added, a UDQ over one of them for instance.
    schedule.insertKeywords(keywordsByStep, parseContext, errors);

    return added;
}

} // namespace Opm
