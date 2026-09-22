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

#ifndef OPM_BLACKOILWELLMODEL_FACILITY_CHECK_IMPL_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_FACILITY_CHECK_IMPL_HEADER_INCLUDED

// Is the state handed to the reservoir linearisation a facility solution at
// this reservoir state? Scored the same way whoever made the decisions:
// the physics (wells, network, limits, nobody throttled without a binding
// limit above) and what one more legacy pass would change. Diagnostic only,
// serial only; the run is left exactly as it was.

#include <fmt/format.h>

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace Opm {

template<typename TypeTag>
void
BlackoilWellModel<TypeTag>::
facilityCheck_(DeferredLogger& deferred_logger)
{
    const int step = simulator_.episodeIndex();
    const auto& schedule = this->schedule();
    const auto& summary_state = this->summaryState();
    const auto& pu = this->phaseUsage();
    const Group& field = schedule.getGroup("FIELD", step);
    const Scalar tol = param_.group_controller_rate_tolerance_;
    const Scalar tol_pressure = param_.group_controller_network_tolerance_;

    auto phase = [&pu](const std::vector<Scalar>& r, const int canonical) {
        return pu.phaseIsActive(canonical) ? r[pu.canonicalToActivePhaseIdx(canonical)] : Scalar(0);
    };
    auto join = [](const std::vector<std::string>& v) {
        std::string s;
        for (const auto& x : v) {
            s += (s.empty() ? "" : ", ") + x;
        }
        return s;
    };

    // The residuals were assembled by the caller at the hand-over state.
    std::vector<std::string> unsolved;
    for (const auto& well : well_container_) {
        if (well->wellIsStopped()) {
            continue;
        }
        if (!well->getWellConvergence(this->groupStateHelper(), B_avg_, /*relax*/ false).converged()) {
            const auto& ws = this->wellState().well(well->indexOfWell());
            unsolved.push_back(well->name() + "(" + (well->isProducer()
                ? WellProducerCMode2String(ws.production_cmode)
                : WellInjectorCMode2String(ws.injection_cmode)) + ")");
        }
    }

    // Own limits of the producers, at the hand-over state.
    std::pair<Scalar, std::string> well_over{Scalar{0}, ""};
    auto worse = [](std::pair<Scalar, std::string>& worst, const Scalar v, const std::string& tag) {
        if (v > worst.first) {
            worst = {v, tag};
        }
    };
    for (const auto& well : well_container_) {
        if (!well->isProducer() || well->wellIsStopped() || !well->wellEcl().predictionMode()) {
            continue;
        }
        const auto& ws = this->wellState().well(well->indexOfWell());
        const auto controls = well->wellEcl().productionControls(summary_state);
        const Scalar oil = -phase(ws.surface_rates, IndexTraits::oilPhaseIdx);
        const Scalar wat = -phase(ws.surface_rates, IndexTraits::waterPhaseIdx);
        const Scalar gas = -phase(ws.surface_rates, IndexTraits::gasPhaseIdx);
        auto rate = [&](const Well::ProducerCMode mode, const Scalar limit, const Scalar current, const char* tag) {
            if (controls.hasControl(mode) && limit > Scalar(0)) {
                worse(well_over, current / limit - Scalar(1), well->name() + ":" + tag);
            }
        };
        rate(Well::ProducerCMode::ORAT, controls.oil_rate, oil, "ORAT");
        rate(Well::ProducerCMode::WRAT, controls.water_rate, wat, "WRAT");
        rate(Well::ProducerCMode::GRAT, controls.gas_rate, gas, "GRAT");
        rate(Well::ProducerCMode::LRAT, controls.liquid_rate, oil + wat, "LRAT");
        if (controls.hasControl(Well::ProducerCMode::BHP) && controls.bhp_limit > Scalar(0)) {
            worse(well_over, Scalar(1) - ws.bhp / controls.bhp_limit, well->name() + ":BHP");
        }
        if (well->wellHasTHPConstraints(summary_state)) {
            const Scalar limit = well->getTHPConstraint(summary_state);
            if (limit > Scalar(0)) {
                worse(well_over, Scalar(1) - ws.thp / limit, well->name() + ":THP");
            }
        }
    }

    // From here the state is written to; everything is put back below.
    const auto saved_active = this->active_wgstate_;
    const auto saved_nupcol = this->nupcol_wgstate_;
    const auto saved_prod_switches = this->switched_prod_groups_;
    const auto saved_inj_switches = this->switched_inj_groups_;
    const auto saved_closed = this->closed_offending_wells_;
    const auto saved_decided = this->controller_decided_wells_;

    std::pair<Scalar, std::string> network_off{Scalar{0}, ""};
    std::pair<Scalar, std::string> group_over{Scalar{0}, ""};
    std::pair<Scalar, std::string> target_move{Scalar{0}, ""};
    std::vector<std::string> idle, group_switches, well_switches, unliftable_switches;
    std::string failure;
    {
        auto log_guard = this->groupStateHelper().pushLogger(/*do_mpi_gather*/ false);
        DeferredLogger scratch;
        try {
            // A legacy pass at the current rates, its switch budgets unspent.
            this->controller_decided_wells_.clear();
            for (const auto& well : well_container_) {
                this->wellState().well(well->indexOfWell()).controller_decided = false;
            }
            this->nupcol_wgstate_ = this->active_wgstate_;
            this->switched_prod_groups_.clear();
            this->switched_inj_groups_.clear();
            this->updateAndCommunicateGroupData(step, /*update_wellgrouptarget*/ true);

            network_off = this->network_.pressureImbalance(step);

            std::set<std::string> binding;
            for (const auto& name : schedule.groupNames(step)) {
                const auto& group = schedule.getGroup(name, step);
                if (!group.isProductionGroup() || !this->groupState().has_production_rates(name)) {
                    continue;
                }
                const auto controls = group.productionControls(summary_state);
                const auto& action = controls.group_limit_action;
                const bool all_rate = action.allRates == Group::ExceedAction::RATE;
                const auto& r = this->groupState().production_rates(name);
                const Scalar oil = phase(r, IndexTraits::oilPhaseIdx);
                const Scalar wat = phase(r, IndexTraits::waterPhaseIdx);
                const Scalar gas = phase(r, IndexTraits::gasPhaseIdx);
                auto limit = [&](const Group::ProductionCMode mode, const Group::ExceedAction act,
                                 const Scalar target, const Scalar current, const char* tag) {
                    if (!group.has_control(mode) || !(all_rate || act == Group::ExceedAction::RATE)) {
                        return;
                    }
                    if (target <= Scalar(0)) {
                        // A zero limit holds everything below it.
                        binding.insert(name);
                        return;
                    }
                    worse(group_over, current / target - Scalar(1), name + ":" + tag);
                    if (current >= (Scalar(1) - tol) * target) {
                        binding.insert(name);
                    }
                };
                limit(Group::ProductionCMode::ORAT, action.oil, controls.oil_target, oil, "ORAT");
                limit(Group::ProductionCMode::WRAT, action.water, controls.water_target, wat, "WRAT");
                limit(Group::ProductionCMode::GRAT, action.gas, controls.gas_target, gas, "GRAT");
                limit(Group::ProductionCMode::LRAT, action.liquid, controls.liquid_target, oil + wat, "LRAT");
            }

            // A well held back by its group needs a binding limit somewhere above it.
            for (const auto& well : well_container_) {
                const auto& ws = saved_active.well_state.well(well->indexOfWell());
                if (!well->isProducer() || well->wellIsStopped()
                    || ws.production_cmode != Well::ProducerCMode::GRUP) {
                    continue;
                }
                bool held = false;
                for (std::string g = well->wellEcl().groupName(); !held; g = schedule.getGroup(g, step).parent()) {
                    held = binding.count(g) > 0;
                    if (g == "FIELD") {
                        break;
                    }
                }
                if (!held) {
                    idle.push_back(well->name());
                }
            }

            std::map<std::string, Group::ProductionCMode> before;
            for (const auto& name : schedule.groupNames(step)) {
                if (this->groupState().has_production_control(name)) {
                    before[name] = this->groupState().production_control(name);
                }
            }
            std::map<std::string, Group::InjectionCMode> inj_before;
            for (const auto& name : schedule.groupNames(step)) {
                for (const Phase ph : {Phase::WATER, Phase::OIL, Phase::GAS}) {
                    if (this->groupState().has_injection_control(name, ph)) {
                        inj_before[name + ":" + std::to_string(static_cast<int>(ph))] = this->groupState().injection_control(name, ph);
                    }
                }
            }
            this->updateGroupControls(field, scratch, step);
            for (const auto& name : schedule.groupNames(step)) {
                for (const Phase ph : {Phase::WATER, Phase::OIL, Phase::GAS}) {
                    if (!this->groupState().has_injection_control(name, ph)) {
                        continue;
                    }
                    const auto now = this->groupState().injection_control(name, ph);
                    const auto it = inj_before.find(name + ":" + std::to_string(static_cast<int>(ph)));
                    if (it == inj_before.end() || it->second != now) {
                        group_switches.push_back(name + ":inj" + std::to_string(static_cast<int>(ph)) + ":"
                            + (it == inj_before.end() ? std::string("-") : Group::InjectionCMode2String(it->second))
                            + "->" + Group::InjectionCMode2String(now));
                    }
                }
            }
            for (const auto& name : schedule.groupNames(step)) {
                if (!this->groupState().has_production_control(name)) {
                    continue;
                }
                const auto now = this->groupState().production_control(name);
                const auto it = before.find(name);
                if (it == before.end() || it->second != now) {
                    group_switches.push_back(name + ":"
                        + (it == before.end() ? std::string("-") : Group::ProductionCMode2String(it->second))
                        + "->" + Group::ProductionCMode2String(now));
                }
            }

            // Legacy's target for the wells it holds on GRUP, against what they produce.
            for (const auto& well : well_container_) {
                const auto& ws = this->wellState().well(well->indexOfWell());
                if (!well->isProducer() || well->wellIsStopped()
                    || ws.production_cmode != Well::ProducerCMode::GRUP || !ws.group_target.has_value()) {
                    continue;
                }
                const auto& r = saved_active.well_state.well(well->indexOfWell()).surface_rates;
                const Scalar oil = -phase(r, IndexTraits::oilPhaseIdx);
                const Scalar wat = -phase(r, IndexTraits::waterPhaseIdx);
                const Scalar gas = -phase(r, IndexTraits::gasPhaseIdx);
                Scalar current = -1;
                switch (ws.group_target->production_cmode) {
                case Group::ProductionCMode::ORAT: current = oil; break;
                case Group::ProductionCMode::WRAT: current = wat; break;
                case Group::ProductionCMode::GRAT: current = gas; break;
                case Group::ProductionCMode::LRAT: current = oil + wat; break;
                default: break;
                }
                const Scalar target = ws.group_target->target_value;
                if (current >= Scalar(0) && std::max(target, current) > Scalar(0)) {
                    worse(target_move, std::abs(target - current) / std::max(target, current), well->name());
                }
            }

            for (const auto& well : well_container_) {
                const auto to = well->legacyWouldSwitch(simulator_, this->groupStateHelper(), this->wellState(),
                                                        tol_pressure);
                if (to.has_value()) {
                    (to->find("(unliftable)") != std::string::npos ? unliftable_switches : well_switches)
                        .push_back(well->name() + ":" + *to);
                }
            }
        } catch (const std::exception& e) {
            failure = e.what();
        }
        log_guard.discard();
    }
    this->active_wgstate_ = saved_active;
    this->nupcol_wgstate_ = saved_nupcol;
    this->switched_prod_groups_ = saved_prod_switches;
    this->switched_inj_groups_ = saved_inj_switches;
    this->closed_offending_wells_ = saved_closed;
    this->controller_decided_wells_ = saved_decided;

    const bool physics_ok = failure.empty() && unsolved.empty() && idle.empty()
        && network_off.first <= tol_pressure && group_over.first <= tol && well_over.first <= tol;
    const bool legacy_ok = failure.empty() && group_switches.empty() && well_switches.empty()
        && target_move.first <= tol && network_off.first <= tol_pressure;

    auto& st = this->facility_check_stats_;
    ++st.checks;
    st.physics_ok += physics_ok;
    st.legacy_ok += legacy_ok;
    st.both_ok += physics_ok && legacy_ok;
    st.physics_only += physics_ok && !legacy_ok;
    st.has_last = true;
    st.last_physics_ok = physics_ok;
    st.last_legacy_ok = legacy_ok;

    const auto& iterCtx = simulator_.problem().iterationContext();
    deferred_logger.debug(fmt::format(
        "Facility check: step {} iteration {}: physics {}, legacy {} | unsolved wells {} [{}] | "
        "network off {:.3f} bar ({}) | group over {:+.1f} % ({}) | well over {:+.1f} % ({}) | "
        "held without a binding limit {} [{}] | legacy would switch groups {} [{}] wells {} [{}], "
        "move a target {:.1f} % ({}); unliftable, not counted {} [{}]{}",
        step, iterCtx.iteration(), physics_ok ? "ok" : "NO", legacy_ok ? "ok" : "NO",
        unsolved.size(), join(unsolved), network_off.first * 1.0e-5, network_off.second,
        100.0 * group_over.first, group_over.second, 100.0 * well_over.first, well_over.second,
        idle.size(), join(idle), group_switches.size(), join(group_switches),
        well_switches.size(), join(well_switches), 100.0 * target_move.first, target_move.second,
        unliftable_switches.size(), join(unliftable_switches),
        failure.empty() ? "" : " | probe failed: " + failure));
}

} // namespace Opm

#endif
