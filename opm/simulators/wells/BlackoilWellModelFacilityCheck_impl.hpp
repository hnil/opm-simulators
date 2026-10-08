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
#include <tuple>
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
        if (well->wellHasTHPConstraints(summary_state)
            && this->controller_off_thp_.count(well->name()) == 0) {
            // WVFPEXP item 4 keeps some wells off their thp on purpose; that is not a violation.
            const Scalar limit = well->getTHPConstraint(summary_state);
            if (limit > Scalar(0)) {
                worse(well_over, Scalar(1) - ws.thp / limit, well->name() + ":THP");
            }
        }
    }

    // Own limits of the injectors: rate, reservoir rate, and bhp and thp from above.
    std::pair<Scalar, std::string> inj_over{Scalar{0}, ""};
    auto injectedPhase = [](const InjectorType t) {
        return t == InjectorType::WATER ? IndexTraits::waterPhaseIdx
             : t == InjectorType::OIL ? IndexTraits::oilPhaseIdx : IndexTraits::gasPhaseIdx;
    };
    for (const auto& well : well_container_) {
        if (!well->isInjector() || well->wellIsStopped() || !well->wellEcl().predictionMode()) {
            continue;
        }
        const auto& ws = this->wellState().well(well->indexOfWell());
        const auto controls = well->wellEcl().injectionControls(summary_state);
        const int ph = injectedPhase(controls.injector_type);
        if (controls.hasControl(Well::InjectorCMode::RATE) && controls.surface_rate > 0.0) {
            worse(inj_over, std::abs(phase(ws.surface_rates, ph)) / controls.surface_rate - Scalar(1),
                  well->name() + ":RATE");
        }
        if (controls.hasControl(Well::InjectorCMode::RESV) && controls.reservoir_rate > 0.0) {
            worse(inj_over, std::abs(phase(ws.reservoir_rates, ph)) / controls.reservoir_rate - Scalar(1),
                  well->name() + ":RESV");
        }
        if (controls.hasControl(Well::InjectorCMode::BHP) && controls.bhp_limit > 0.0) {
            worse(inj_over, ws.bhp / controls.bhp_limit - Scalar(1), well->name() + ":BHP");
        }
        if (well->wellHasTHPConstraints(summary_state)) {
            const Scalar limit = well->getTHPConstraint(summary_state);
            if (limit > Scalar(0)) {
                worse(inj_over, ws.thp / limit - Scalar(1), well->name() + ":THP");
            }
        }
    }

    // Injection groups, per phase, with legacy's own sums: RATE, RESV, REIN (of the reinjection
    // group's production), VREP (of the voidage group's reservoir voidage).
    std::pair<Scalar, std::string> group_inj_over{Scalar{0}, ""};
    {
        const auto& helper = this->groupStateHelper();
        auto pos = [&pu](const int canonical) {
            return pu.phaseIsActive(canonical) ? pu.canonicalToActivePhaseIdx(canonical) : -1;
        };
        auto resSum = [&](const Group& g, const bool injector) {
            Scalar s = 0;
            for (const int c : {IndexTraits::waterPhaseIdx, IndexTraits::oilPhaseIdx, IndexTraits::gasPhaseIdx}) {
                if (pos(c) >= 0) { s += helper.sumWellResRates(g, pos(c), injector); }
            }
            return this->comm().sum(s);
        };
        for (const auto& name : schedule.groupNames(step)) {
            const auto& group = schedule.getGroup(name, step);
            if (!group.isInjectionGroup()) {
                continue;
            }
            for (const auto& [ph, canonical, tag] : {std::tuple{Phase::WATER, IndexTraits::waterPhaseIdx, "WAT"},
                                                     std::tuple{Phase::GAS, IndexTraits::gasPhaseIdx, "GAS"},
                                                     std::tuple{Phase::OIL, IndexTraits::oilPhaseIdx, "OIL"}}) {
                if (!group.hasInjectionControl(ph) || pos(canonical) < 0) {
                    continue;
                }
                const auto controls = group.injectionControls(ph, summary_state);
                const Scalar injected = this->comm().sum(helper.sumWellSurfaceRates(group, pos(canonical), true));
                auto over = [&](const Scalar current, const Scalar target, const char* mode) {
                    if (target > Scalar(0)) {
                        worse(group_inj_over, current / target - Scalar(1), name + ":" + tag + ":" + mode);
                    }
                };
                if (group.has_control(ph, Group::InjectionCMode::RATE)) {
                    over(injected, group.has_gpmaint_control(ph, Group::InjectionCMode::RATE)
                                   ? this->groupState().gpmaint_target(name) : controls.surface_max_rate, "RATE");
                }
                if (group.has_control(ph, Group::InjectionCMode::RESV)) {
                    over(this->comm().sum(helper.sumWellResRates(group, pos(canonical), true)),
                         group.has_gpmaint_control(ph, Group::InjectionCMode::RESV)
                             ? this->groupState().gpmaint_target(name) : controls.resv_max_rate, "RESV");
                }
                if (group.has_control(ph, Group::InjectionCMode::REIN)) {
                    const Group& from = schedule.getGroup(controls.reinj_group, step);
                    const Scalar produced = this->comm().sum(helper.sumWellSurfaceRates(from, pos(canonical), false));
                    over(injected, controls.target_reinj_fraction * produced, "REIN");
                }
                if (group.has_control(ph, Group::InjectionCMode::VREP)) {
                    const Group& from = schedule.getGroup(controls.voidage_group, step);
                    over(resSum(group, true), controls.target_void_fraction * resSum(from, false), "VREP");
                }
            }
        }
    }

    // The route's answer against what each well then did; exact hand-over gives zero.
    std::pair<Scalar, std::string> off_assigned{Scalar{0}, ""};
    for (const auto& [name, a] : this->controller_assigned_rates_) {
        if (!this->wellState().has(name)) {
            continue;
        }
        const auto& q = this->wellState().well(name).surface_rates;
        for (const auto& [c, tag, floor] : {std::tuple{IndexTraits::oilPhaseIdx, "oil", 1.0},
                                            std::tuple{IndexTraits::waterPhaseIdx, "wat", 1.0},
                                            std::tuple{IndexTraits::gasPhaseIdx, "gas", 1000.0}}) {
            const Scalar x = std::abs(phase(a, c)), y = std::abs(phase(q, c));
            if (std::max(x, y) > Scalar(floor / 86400.0)) {
                worse(off_assigned, std::abs(x - y) / std::max(x, y), name + ":" + tag);
            }
        }
    }

    // Wells the run closed while the schedule has them open: no well object left to probe.
    std::vector<std::string> closed_by_run;
    for (const auto& name : schedule.wellNames(step)) {
        const auto& w = schedule.getWell(name, step);
        if (w.isProducer() && w.getStatus() == Well::Status::OPEN && w.predictionMode()
            && this->wellState().has(name) && this->wellState().well(name).status == Well::Status::SHUT) {
            closed_by_run.push_back(name);
        }
    }

    // From here the state is written to; everything is put back below.
    const auto saved_active = this->active_wgstate_;
    const auto saved_nupcol = this->nupcol_wgstate_;
    const auto saved_prod_switches = this->switched_prod_groups_;
    const auto saved_inj_switches = this->switched_inj_groups_;
    const auto saved_closed = this->closed_offending_wells_;
    const auto saved_decided = this->controller_decided_wells_;
    const auto saved_inj_decided = this->controller_injection_decided_;

    std::pair<Scalar, std::string> network_off{Scalar{0}, ""};
    std::pair<Scalar, std::string> group_over{Scalar{0}, ""};
    std::pair<Scalar, std::string> target_move{Scalar{0}, ""};
    std::vector<std::string> idle, group_switches, well_switches, unliftable_switches, liftable, closed_would_flow,
        kept_would_flow, cliffs, unconfirmed;
    std::string failure;
    {
        auto log_guard = this->groupStateHelper().pushLogger(/*do_mpi_gather*/ false);
        DeferredLogger scratch;
        // Each stage runs everywhere or nowhere: a rank that throws before a collective would
        // otherwise leave the others waiting in it.
        auto stage = [&](auto&& body) {
            int failed = failure.empty() ? 0 : 1;
            if (failed == 0) {
                try {
                    body();
                } catch (const std::exception& e) {
                    failure = e.what();
                    failed = 1;
                }
            }
            return this->comm().max(failed) == 0;
        };
        stage([&] {
            // A legacy pass at the current rates, its switch budgets unspent.
            this->controller_decided_wells_.clear();
            this->controller_injection_decided_.clear();
            for (const auto& well : well_container_) {
                this->wellState().well(well->indexOfWell()).controller_decided = false;
            }
            this->nupcol_wgstate_ = this->active_wgstate_;
            this->switched_prod_groups_.clear();
            this->switched_inj_groups_.clear();
        })
        && stage([&] {
            this->updateAndCommunicateGroupData(step, /*update_wellgrouptarget*/ true);

            network_off = this->network_.pressureImbalance(step);
            if (std::getenv("OPM_FACILITY_CHECK_DETAIL") && network_off.first > 0.1 * unit::barsa
                && schedule.hasGroup(network_off.second, step)) {
                std::string rows;
                for (const auto& wname : schedule.getGroup(network_off.second, step).wells()) {
                    if (!this->wellState().has(wname)) { continue; }
                    const auto& ws = this->wellState().well(wname);
                    const auto it = this->controller_assigned_rates_.find(wname);
                    auto ph = [&](const std::vector<Scalar>& r, const int c) {
                        return pu.phaseIsActive(c) ? -r[pu.canonicalToActivePhaseIdx(c)] * 86400.0 : 0.0;
                    };
                    rows += fmt::format(" | {} cmode {} status {} thp {:.2f} bhp {:.2f} alq {:.0f}: oil {:.1f} gas {:.0f} wat {:.1f}", wname,
                                        static_cast<int>(ws.production_cmode), static_cast<int>(ws.status),
                                        ws.thp / unit::barsa, ws.bhp / unit::barsa, ws.alq_state.get() * 86400.0,
                                        ph(ws.surface_rates, IndexTraits::oilPhaseIdx), ph(ws.surface_rates, IndexTraits::gasPhaseIdx),
                                        ph(ws.surface_rates, IndexTraits::waterPhaseIdx));
                    auto tgt = [](const auto& t) {
                        return t.has_value() ? fmt::format("{} {:.1f}", Group::ProductionCMode2String(t->production_cmode),
                                                           t->target_value * 86400.0)
                                             : std::string("none");
                    };
                    rows += fmt::format(" target {} fallback {}{}{}", tgt(ws.group_target), tgt(ws.group_target_fallback),
                                        ws.use_group_target_fallback ? " (using fallback)" : "",
                                        ws.controller_decided ? " decided" : "");
                    if (it != this->controller_assigned_rates_.end()) {
                        rows += fmt::format(" (assigned oil {:.1f} gas {:.0f} wat {:.1f})", ph(it->second, IndexTraits::oilPhaseIdx),
                                            ph(it->second, IndexTraits::gasPhaseIdx), ph(it->second, IndexTraits::waterPhaseIdx));
                    }
                }
                // What the check's network computation is fed for this node: the group's rates, and its lift gas.
                std::string group_rates;
                if (this->groupState().has_production_rates(network_off.second)) {
                    const auto& g = this->groupState().production_rates(network_off.second);
                    auto at = [&](const int c) { return pu.phaseIsActive(c) ? g[pu.canonicalToActivePhaseIdx(c)] * 86400.0 : 0.0; };
                    group_rates = fmt::format(" | group rates oil {:.1f} gas {:.0f} wat {:.1f}", at(IndexTraits::oilPhaseIdx),
                                              at(IndexTraits::gasPhaseIdx), at(IndexTraits::waterPhaseIdx));
                }
                double lift = 0.0;
                for (const auto& wname : schedule.getGroup(network_off.second, step).wells()) {
                    if (this->wellState().has(wname) && this->wellState().isOpen(wname)) {
                        lift += this->wellState().well(wname).alq_state.get() * 86400.0;
                    }
                }
                const auto& stored = this->network_.nodePressures();
                const auto sp = stored.find(network_off.second);
                OpmLog::debug(fmt::format("Facility check detail: {} off {:.3f} bar, stored {:.3f} bar{} lift {:.0f}{}",
                                          network_off.second, network_off.first / unit::barsa,
                                          sp != stored.end() ? sp->second / unit::barsa : -1.0, group_rates, lift, rows));
            }
        })
        && stage([&] {

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
        })
        && stage([&] {
            // Zero rate is a solution only where the well cannot flow at its thp limit: the stable crossing of
            // its own inflow and tubing, on a well made for it, without a solve. A well the run closed may have
            // been closed on an economic limit, so it is listed apart.
            for (const auto& name : schedule.wellNames(step)) {
                const auto& w = schedule.getWell(name, step);
                if (!w.isProducer() || w.getStatus() != Well::Status::OPEN || !w.predictionMode()
                    || !this->wellState().has(name)) {
                    continue;
                }
                const auto live = std::find_if(well_container_.begin(), well_container_.end(),
                                               [&name](const auto& wp) { return wp->name() == name; });
                const bool stopped = live != well_container_.end() && (*live)->wellIsStopped();
                const bool closed = live == well_container_.end()
                    && this->wellState().well(name).status == Well::Status::SHUT;
                if (!stopped && !closed) {
                    continue;
                }
                auto probe = this->createWellForWellTest(name, step, scratch);
                probe->init(depth_, gravity_, B_avg_, true);
                probe->setWellEfficiencyFactor(stopped ? (*live)->wellEfficiencyFactor() : Scalar(1));
                probe->setVFPProperties(this->vfp_properties_.get());
                probe->setGuideRate(&this->guideRate_);
                if (probe->isVFPActive(scratch)) {
                    probe->setPrevSurfaceRates(this->wellState(), this->prevWellState());
                }
                this->network_.initializeWell(*probe);
                probe->setStableThpCrossing(true);
                probe->calculateExplicitQuantities(simulator_, this->groupStateHelper());
                const Scalar alq = probe->getALQ(this->wellState());
                std::optional<Scalar> bhp;
                auto flowsAt = [&](const std::optional<Scalar> thp) -> std::optional<std::vector<Scalar>> {
                    if (thp.has_value()) {
                        probe->setDynamicThpLimit(*thp);
                        bhp = probe->computeBhpAtThpLimitProdWithAlq(simulator_, this->groupStateHelper(), summary_state,
                                                                     alq, /*iterate_if_no_solution*/ false);
                    } else {
                        bhp = WellBhpThpCalculator<Scalar, IndexTraits>(*probe).mostStrictBhpFromBhpLimits(summary_state);
                    }
                    if (!bhp.has_value()) {
                        return std::nullopt;
                    }
                    std::vector<Scalar> flux(probe->numConservationQuantities(), Scalar(0));
                    probe->computeWellRatesWithBhp(simulator_, *bhp, flux, scratch);
                    // In active phase order already.
                    const auto& surface = flux;
                    if (std::none_of(surface.begin(), surface.end(), [](const Scalar v) { return v < Scalar(0); })) {
                        return std::nullopt;
                    }
                    return surface;
                };
                const bool has_thp = probe->wellHasTHPConstraints(summary_state);
                const auto q = flowsAt(has_thp ? std::optional<Scalar>(probe->getTHPConstraint(summary_state))
                                               : std::nullopt);
                // The crossing is found on the inflow without a solve; the well's own equations at that bhp
                // must agree, or it is not counted.
                // Stopped before this step and no WTEST: the deck keeps it stopped; reported, not counted.
                const bool earlier = stopped && this->prevWellState().has(name)
                    && this->prevWellState().well(name).status != Well::Status::OPEN;
                Scalar full_model = 0;
                if (q.has_value() && has_thp && stopped && !earlier) {
                    full_model = probe->thpMarginWithIterations(simulator_, this->groupStateHelper(), *bhp);
                    if (full_model < -tol_pressure) {
                        unconfirmed.push_back(fmt::format("{}@{:.1f} oil {:.0f}: {:+.2f} bar", name,
                                                          probe->getTHPConstraint(summary_state) * 1.0e-5,
                                                          -phase(*q, IndexTraits::oilPhaseIdx) * 86400.0,
                                                          full_model * 1.0e-5));
                        continue;
                    }
                }
                if (q.has_value()) {
                    const Scalar thp = has_thp ? probe->getTHPConstraint(summary_state) : Scalar(0);
                    // A cliff when it cannot flow either at the node pressure its own flow gives: neither
                    // state is a solution.
                    std::string at_open;
                    bool cliff = false;
                    if (stopped && !earlier && has_thp) {
                        auto& ws = this->wellState().well(name);
                        const auto saved_ws = ws;
                        ws.open();
                        ws.surface_rates = *q;
                        this->updateAndCommunicateGroupData(step, /*update_wellgrouptarget*/ false);
                        const auto p_open = this->network_.pressuresAtGroupRates(step);
                        ws = saved_ws;
                        this->updateAndCommunicateGroupData(step, /*update_wellgrouptarget*/ false);
                        if (const auto it = p_open.find(w.groupName()); it != p_open.end() && it->second > thp) {
                            const auto q_open = flowsAt(it->second);
                            cliff = !q_open.has_value();
                            at_open = fmt::format(" (open {:.1f}: {})", it->second * 1.0e-5,
                                                  cliff ? std::string("dead") : fmt::format("oil {:.0f}",
                                                      -phase(*q_open, IndexTraits::oilPhaseIdx) * 86400.0));
                        }
                    }
                    (!stopped ? closed_would_flow : earlier ? kept_would_flow : cliff ? cliffs : liftable).push_back(
                        fmt::format("{}@{:.1f} oil {:.0f} gas {:.0f}{}", name, thp * 1.0e-5,
                                    -phase(*q, IndexTraits::oilPhaseIdx) * 86400.0,
                                    -phase(*q, IndexTraits::gasPhaseIdx) * 86400.0, at_open));
                }
            }
        });
        log_guard.discard();
    }
    this->active_wgstate_ = saved_active;
    this->nupcol_wgstate_ = saved_nupcol;
    this->switched_prod_groups_ = saved_prod_switches;
    this->switched_inj_groups_ = saved_inj_switches;
    this->closed_offending_wells_ = saved_closed;
    this->controller_decided_wells_ = saved_decided;
    this->controller_injection_decided_ = saved_inj_decided;

    if (this->comm().size() > 1) {
        // Each rank sees its own wells: the check is the worst of them, and the counts are the totals.
        for (auto* p : {&network_off, &group_over, &target_move, &well_over, &inj_over, &group_inj_over, &off_assigned}) {
            p->first = this->comm().max(p->first);
        }
        const int any_failure = this->comm().max(failure.empty() ? 0 : 1);
        if (any_failure != 0 && failure.empty()) {
            failure = "another rank";
        }
        const auto counts = [this](std::vector<std::string>& v) {
            const int n = this->comm().sum(static_cast<int>(v.size()));
            if (n > 0 && v.empty()) { v.emplace_back("(another rank)"); }
        };
        for (auto* v : {&unsolved, &idle, &group_switches, &well_switches, &unliftable_switches, &liftable,
                        &closed_by_run, &closed_would_flow, &kept_would_flow, &cliffs,
                        &unconfirmed}) {
            counts(*v);
        }
    }

    const bool physics_ok = failure.empty() && unsolved.empty() && idle.empty()
        && network_off.first <= tol_pressure && group_over.first <= tol && well_over.first <= tol
        && inj_over.first <= tol && group_inj_over.first <= tol;
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
        "injector over {:+.1f} % ({}) | injection group over {:+.1f} % ({}) | "
        "held without a binding limit {} [{}] | legacy would switch groups {} [{}] wells {} [{}], "
        "move a target {:.1f} % ({}); unliftable, not counted {} [{}] | stopped this step but liftable {} [{}], "
        "on a cliff {} [{}], kept stopped, would flow {} [{}], not confirmed by the well's equations {} [{}] | "
        "closed by the run {} [{}], would flow {} [{}] | off assignment {:.1f} % ({}){}",
        step, iterCtx.iteration(), physics_ok ? "ok" : "NO", legacy_ok ? "ok" : "NO",
        unsolved.size(), join(unsolved), network_off.first * 1.0e-5, network_off.second,
        100.0 * group_over.first, group_over.second, 100.0 * well_over.first, well_over.second,
        100.0 * inj_over.first, inj_over.second, 100.0 * group_inj_over.first, group_inj_over.second,
        idle.size(), join(idle), group_switches.size(), join(group_switches),
        well_switches.size(), join(well_switches), 100.0 * target_move.first, target_move.second,
        unliftable_switches.size(), join(unliftable_switches), liftable.size(), join(liftable),
        cliffs.size(), join(cliffs), kept_would_flow.size(), join(kept_would_flow), unconfirmed.size(), join(unconfirmed),
        closed_by_run.size(), join(closed_by_run), closed_would_flow.size(), join(closed_would_flow),
        100.0 * off_assigned.first, off_assigned.second,
        failure.empty() ? "" : " | probe failed: " + failure));
}

} // namespace Opm

#endif
