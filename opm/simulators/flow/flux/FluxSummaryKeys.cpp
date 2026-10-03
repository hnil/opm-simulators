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

#include <opm/simulators/flow/flux/FluxSummaryKeys.hpp>

#include <opm/input/eclipse/Schedule/Action/ActionX.hpp>
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/RequisiteSummaryVector.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>

#include <opm/io/eclipse/SummaryNode.hpp>

#include <algorithm>
#include <initializer_list>
#include <unordered_set>

namespace Opm {

std::vector<std::string> fluxSummaryKeywords(const Schedule& schedule,
                                             const bool oil,
                                             const bool water,
                                             const bool gas)
{
    auto keywords = std::unordered_set<std::string>{};

    const auto addAll = [&keywords](std::initializer_list<const char*> names)
    {
        for (const auto* name : names) {
            keywords.insert(name);
        }
    };

    // Surface rates and cumulatives for the conserved quantities of every
    // active phase. A USEFLUX run needs these to reconstruct the contribution
    // of wells that fall outside the sector.
    if (oil) {
        addAll({"WOPR", "WOPT", "WOIR", "WOIT"});
    }
    if (water) {
        addAll({"WWPR", "WWPT", "WWIR", "WWIT"});
    }
    if (gas) {
        addAll({"WGPR", "WGPT", "WGIR", "WGIT"});
    }

    // Reservoir volume rates and cumulatives.
    addAll({"WVPR", "WVPT", "WVIR", "WVIT"});

    // Everything referenced by the deck's UDQ DEFINE expressions and by ACTIONX
    // conditions, so that both can be evaluated through the standard code paths
    // in the reduced run.
    for (const auto& udq : schedule.unique<UDQConfig>()) {
        udq.second.required_summary(keywords);

        // A definition's requirements stop at the summary vectors it reads;
        // the UDQs it reads are filtered out, on the grounds that the run
        // computes those. It can only compute one whose inputs it still has,
        // though, and an ASSIGN has no inputs at all, so carry them too.
        udq.second.requiredUDQs(keywords);
    }

    for (const auto& action : schedule.back().actions.get()) {
        action.required_summary(keywords);
    }

    auto sorted = std::vector<std::string>(keywords.begin(), keywords.end());
    sorted.erase(std::remove_if(sorted.begin(), sorted.end(),
                                [](const std::string& keyword)
                                { return keyword.empty(); }),
                 sorted.end());

    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());

    return sorted;
}

std::string_view fluxSummaryKeywordOf(std::string_view key)
{
    return key.substr(0, key.find(':'));
}

std::vector<std::string> fluxSummaryKeys(const Schedule& schedule,
                                         const bool oil,
                                         const bool water,
                                         const bool gas,
                                         const std::vector<std::string>& availableKeys)
{
    const auto keywords = fluxSummaryKeywords(schedule, oil, water, gas);
    const auto wanted = std::unordered_set<std::string>(keywords.begin(), keywords.end());

    // The keywords are bare, so expand the well and group level ones over the
    // objects they can apply to.
    const auto& wells = schedule.wellNames();
    const auto& groups = schedule.groupNames();

    auto keys = std::vector<std::string>{};
    keys.reserve(keywords.size());

    for (const auto& keyword : keywords) {
        switch (EclIO::SummaryNode::category_from_keyword(keyword)) {
        case EclIO::SummaryNode::Category::Well:
            for (const auto& well : wells) {
                keys.push_back(keyword + ':' + well);
            }
            break;

        case EclIO::SummaryNode::Category::Group:
            for (const auto& group : groups) {
                keys.push_back(keyword + ':' + group);
            }
            break;

        case EclIO::SummaryNode::Category::Field:
        case EclIO::SummaryNode::Category::Miscellaneous:
            keys.push_back(keyword);
            break;

        default:
            // Region, block, connection, segment, aquifer and node level
            // quantities name an object the schedule cannot enumerate: nothing
            // in it says which regions of which region set the run reports on.
            // They are picked out of availableKeys below instead.
            break;
        }
    }

    // Whatever the run itself holds under a wanted keyword. This is the only
    // way the categories above can be reached, and leaving them out is not the
    // harmless omission it looks: a reduced run evaluating an ACTIONX that
    // names RPR__REC:3 aborts with
    //
    //   Summary vector RPR__REC:3 is unknown
    //
    // and no indication that the vector should have come from the parent.
    for (const auto& key : availableKeys) {
        switch (EclIO::SummaryNode::category_from_keyword(key)) {
        case EclIO::SummaryNode::Category::Well:
        case EclIO::SummaryNode::Category::Group:
        case EclIO::SummaryNode::Category::Field:
        case EclIO::SummaryNode::Category::Miscellaneous:
            // Already covered by the expansion above, which does not depend on
            // the run having got far enough to hold a value yet.
            continue;

        default:
            break;
        }

        if (wanted.count(std::string{fluxSummaryKeywordOf(key)}) != 0) {
            keys.push_back(key);
        }
    }

    // And whatever the deck's expressions name outright. A condition on
    // RPR__REC:3 is answered by the parent whether or not its SUMMARY section
    // ever mentioned region 3, so the key can be absent from availableKeys and
    // still be needed.
    auto named = RequisiteSummaryVectors{};
    for (const auto& udq : schedule.unique<UDQConfig>()) {
        udq.second.requisiteSummaryVectors(named);
    }
    for (const auto& action : schedule.back().actions.get()) {
        action.requisiteSummaryVectors(named);
    }

    for (const auto& vector : named) {
        switch (EclIO::SummaryNode::category_from_keyword(vector.keyword)) {
        case EclIO::SummaryNode::Category::Region:
        case EclIO::SummaryNode::Category::Segment:
        case EclIO::SummaryNode::Category::Node:
            break;

        default:
            // Wells and groups are expanded above. A block or connection names
            // its cell by I, J and K while the summary names it by global
            // index, so the two spellings cannot be reconciled here.
            continue;
        }

        auto key = vector.keyword;
        for (const auto& argument : vector.arguments) {
            key += ':';
            key += argument;
        }

        keys.push_back(std::move(key));
    }

    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());

    return keys;
}

} // namespace Opm
