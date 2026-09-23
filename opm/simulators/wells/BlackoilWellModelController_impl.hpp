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

#ifndef OPM_BLACKOILWELLMODEL_CONTROLLER_IMPL_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_CONTROLLER_IMPL_HEADER_INCLUDED

// The group controller's decide step with a production network: every
// well's IPR, the node pressures and the group tree in one system, solved
// by the reduced route, scored by the judge, and handed to the wells as
// controls and targets. Without a network the balancer decides instead.

#include <opm/input/eclipse/Schedule/Group/GSatProd.hpp>
#include <opm/input/eclipse/Schedule/Network/Branch.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/Network/Node.hpp>
#include <opm/input/eclipse/Schedule/ScheduleState.hpp>
#include <opm/input/eclipse/Schedule/VFPProdTable.hpp>

#include <opm/simulators/wells/ProdGroupTreeBalancer.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/network/NetworkJudge.hpp>
#include <opm/simulators/wells/network/NetworkProductionSystem.hpp>
#include <opm/simulators/wells/network/NetworkReducedSolve.hpp>
#include <opm/simulators/wells/network/NetworkTubingExtension.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <array>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace Opm {

template<typename TypeTag>
void
BlackoilWellModel<TypeTag>::
controllerDecide_(DeferredLogger& deferred_logger, const bool write_rates)
{
    if (param_.enable_group_controller_network_ && controllerNetworkDecide_(deferred_logger)) {
        this->controller_network_owned_ = true;
    } else {
        this->controller_network_owned_ = false;
        runControllerBalance_(prepareWellsForBalancing_(deferred_logger), deferred_logger, write_rates);
    }
}

template<typename TypeTag>
std::string
BlackoilWellModel<TypeTag>::
controllerDecisionSignature_() const
{
    std::string sig;
    for (const auto& [name, cmode] : this->controller_assigned_cmode_) {
        Scalar q = 0;
        if (const auto it = this->controller_assigned_rates_.find(name); it != this->controller_assigned_rates_.end()) {
            for (const auto v : it->second) { q += std::abs(v); }
        }
        // Rates in bins of the rate tolerance, pressures in steps of the network tolerance.
        const long bin = q > Scalar{0}
            ? std::lround(std::log(q * 86400.0 + 1.0) / std::log1p(param_.group_controller_rate_tolerance_)) : 0;
        sig += fmt::format("{}:{}:{};", name, static_cast<int>(cmode), bin);
    }
    for (const auto& [node, p] : this->network_.nodePressures()) {
        sig += fmt::format("{}={};", node, std::lround(p / param_.group_controller_network_tolerance_));
    }
    return sig;
}

template<typename TypeTag>
void
BlackoilWellModel<TypeTag>::
controllerMarkDecided_()
{
    for (const auto& well : well_container_) {
        this->wellState().well(well->indexOfWell()).controller_decided =
            this->controller_decided_wells_.count(well->name()) > 0
            || this->controller_injection_decided_.count(well->name()) > 0;
    }
}

template<typename TypeTag>
void
BlackoilWellModel<TypeTag>::
controllerRefreshIpr_(DeferredLogger& deferred_logger)
{
    // The system needs every well's rate response to its own bhp: the implicit
    // IPR, which the well solve keeps only where its own logic needs it. At zero
    // rate that linearisation can be singular; only a degenerate one -- steeper
    // than the connections allow -- is replaced, by the IPR from the last time
    // the well flowed, or the connections' where it never has.
    auto& cache = this->controller_last_flowing_ipr_;
    for (const auto& well : well_container_) {
        if (!well->wellEcl().predictionMode() || !well->isProducer()) {
            continue;
        }
        auto& ws = this->wellState().well(well->indexOfWell());
        if (ws.status != WellStatus::OPEN) {
            continue;
        }
        well->updateIPRImplicit(simulator_, this->groupStateHelper(), this->wellState());
        const bool zero = std::all_of(ws.surface_rates.begin(), ws.surface_rates.end(),
                                      [](const Scalar q) { return q == Scalar{0}; });
        // The implicit intercept is the well's unconverged rate, taken before the reservoir
        // moved (5 bar of shut-in pressure on STDW's C-2H); the true inflow at its bhp is current.
        static const bool true_intercept = [] {
            const char* v = std::getenv("OPM_CONTROLLER_IPR_INTERCEPT");
            return v == nullptr || std::string(v) != "stale";
        }();
        if (!zero && true_intercept) {
            std::vector<Scalar> q(ws.implicit_ipr_a.size(), Scalar{0});
            well->computeWellRatesWithBhp(simulator_, ws.bhp, q, deferred_logger);
            for (std::size_t ph = 0; ph < q.size(); ++ph) {
                ws.implicit_ipr_a[ph] = ws.implicit_ipr_b[ph] * ws.bhp - q[ph];
            }
        }
        auto& source = this->controller_ipr_source_[well->name()];
        if (!zero) {
            cache[well->name()] = {ws.implicit_ipr_a, ws.implicit_ipr_b};
            source = "flowing";
            continue;
        }
        const auto implicit_a = ws.implicit_ipr_a, implicit_b = ws.implicit_ipr_b;
        well->setImplicitIprFromConnections(simulator_, this->wellState(), deferred_logger);
        auto steepest = [](const std::vector<Scalar>& b) {
            Scalar m = 0;
            for (const auto v : b) { m = std::max(m, std::abs(v)); }
            return m;
        };
        const bool degenerate = steepest(implicit_b) > Scalar{10} * steepest(ws.implicit_ipr_b);
        // OPM_CONTROLLER_REVIVAL_IPR=zero-rate keeps the old choice: the tangent at zero rate
        // unless it is degenerate. It can be several times flatter than the flowing well's,
        // and a well the route shut on one Newton iterate then never comes back.
        static const bool last_flowing = [] {
            const char* v = std::getenv("OPM_CONTROLLER_REVIVAL_IPR");
            return v == nullptr || std::string(v) != "zero-rate";
        }();
        const auto it = cache.find(well->name());
        if (last_flowing && it != cache.end()) {
            // The productivity it had when it last flowed, at the static pressure it sees
            // now: each phase's line keeps its slope and is moved to today's zero-inflow bhp.
            const auto& now_a = degenerate ? ws.implicit_ipr_a : implicit_a;
            const auto& now_b = degenerate ? ws.implicit_ipr_b : implicit_b;
            const int oil = this->phaseUsage().canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
            const Scalar p_static = (oil >= 0 && now_b[oil] > Scalar{0}) ? now_a[oil] / now_b[oil] : Scalar{-1};
            ws.implicit_ipr_a = it->second.first;
            ws.implicit_ipr_b = it->second.second;
            if (p_static > Scalar{0}) {
                for (std::size_t ph = 0; ph < ws.implicit_ipr_b.size(); ++ph) {
                    ws.implicit_ipr_a[ph] = ws.implicit_ipr_b[ph] * p_static;
                }
            }
            source = "last flowing slope at today's static pressure";
        } else if (!degenerate) {
            ws.implicit_ipr_a = implicit_a;
            ws.implicit_ipr_b = implicit_b;
            source = "zero-rate tangent";
        } else if (it != cache.end()) {
            ws.implicit_ipr_a = it->second.first;
            ws.implicit_ipr_b = it->second.second;
            source = "last flowing";
        } else {
            source = "connections";
        }
    }
}

template<typename TypeTag>
bool
BlackoilWellModel<TypeTag>::
controllerThpRouteApplies_() const
{
    if (!param_.enable_group_controller_thp_route_ || !param_.enable_group_controller_network_) {
        return false;
    }
    const int step = simulator_.episodeIndex();
    for (const auto& name : this->schedule().wellNames(step)) {
        const auto& well = this->schedule().getWell(name, step);
        if (!well.isProducer() || !well.predictionMode() || well.getStatus() == Well::Status::SHUT) {
            continue;
        }
        const auto controls = well.productionControls(this->summaryState());
        if (controls.vfp_table_number > 0 && controls.thp_limit > 0.0) {
            return true;
        }
    }
    return false;
}

template<typename TypeTag>
bool
BlackoilWellModel<TypeTag>::
controllerNetworkDecide_(DeferredLogger& deferred_logger)
{
    this->controller_rejected_ = false;
    OPM_TIMEFUNCTION();
    using Sys = NetworkSolve::ProductionSystem<Scalar>;
    using Ctrl = typename Sys::Control;
    using Mode = typename Sys::Mode;
    const int reportStepIdx = simulator_.episodeIndex();
    const auto& schedule = this->schedule();
    const auto& network = schedule[reportStepIdx].network();
    // Without a network the route still carries the wells that have a tubing table: each
    // with its own thp limit where a network would give its node's pressure. The capacity
    // the tree is handed and the control the well is given then come from one model.
    const bool no_network = !network.active();
    if (no_network && !controllerThpRouteApplies_()) {
        return false;
    }
    if (this->comm().size() > 1) {
        OPM_DEFLOG_THROW(std::runtime_error,
                         "The group controller's network route is serial only for now", deferred_logger);
    }
    auto giveUp = [&](const std::string& why) {
        deferred_logger.debug(fmt::format("Controller: the network route is not taken at report step {} ({}); "
                                          "the balancer decides and legacy balances the network", reportStepIdx, why));
        return false;
    };
    const auto& units = schedule.getUnits();
    const auto& summary_state = this->summaryState();
    const auto& pu = this->phaseUsage();
    // Water, oil, gas -- the order VFPPROD wants -- as positions in the active-phase arrays.
    const std::array<int, Sys::NP> pos{
        pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx),
        pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx),
        pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)};
    if (std::any_of(pos.begin(), pos.end(), [](const int p) { return p < 0; })) {
        return giveUp("the network needs all three phases");
    }
    controllerRefreshIpr_(deferred_logger);

    // What the deck allows, not what legacy has the well on: available for group
    // control with a rate target somewhere above it.
    const auto groupControllable = [&](const Well& well) {
        if (!well.isAvailableForGroupControl()) {
            return false;
        }
        for (std::string g = well.groupName(); !g.empty(); ) {
            const auto& grp = schedule.getGroup(g, reportStepIdx);
            if (grp.isProductionGroup()) {
                const auto ctl = grp.productionControls(summary_state);
                using C = Group::ProductionCMode;
                const bool rate = ctl.cmode == C::ORAT || ctl.cmode == C::LRAT
                               || ctl.cmode == C::GRAT || ctl.cmode == C::WRAT;
                if (rate && grp.has_control(ctl.cmode)) {
                    return true;
                }
            }
            if (grp.parent() == g) {
                break;
            }
            g = grp.parent();
        }
        return false;
    };
    // Diagnostic: OPM_CONTROLLER_IMPLICIT_GUIDES reads the guide rates afresh at every decision.
    static const bool explicit_guides = std::getenv("OPM_CONTROLLER_IMPLICIT_GUIDES") == nullptr;
    if (const std::pair<double, double> key{simulator_.time(), simulator_.timeStepSize()};
        key != this->controller_step_guides_key_) {
        this->controller_step_guides_key_ = key;
        this->controller_step_guides_.clear();
    }
    // The deck's guide rate on the mode of the nearest targeted group above the well.
    // The deck's guide rate of a well on a given phase, the time step's value.
    const auto deckGuideOn = [&](const std::string& name, const GuideRateModel::Target target) -> Scalar {
        const std::string key = name + ":" + std::to_string(static_cast<int>(target));
        if (explicit_guides) {
            if (const auto it = this->controller_step_guides_.find(key); it != this->controller_step_guides_.end()) {
                return it->second;
            }
        }
        const auto& gr = this->guideRate();
        if (!gr.has(name) && !gr.hasPotentials(name)) {
            return Scalar{0};
        }
        const auto value = static_cast<Scalar>(gr.get(name, target, this->groupStateHelper().getWellRateVector(name)));
        if (explicit_guides && value > Scalar{0}) {
            this->controller_step_guides_[key] = value;
        }
        return value;
    };
    const auto targetOf = [](const Mode m) {
        switch (m) {
        case Mode::Water:  return GuideRateModel::Target::WAT;
        case Mode::Gas:    return GuideRateModel::Target::GAS;
        case Mode::Liquid: return GuideRateModel::Target::LIQ;
        case Mode::Resv:   return GuideRateModel::Target::RES;
        default:           return GuideRateModel::Target::OIL;
        }
    };
    // First guess before the tree exists: the mode of the nearest targeted group above the well.
    const auto deckGuide = [&](const std::string& name) -> Scalar {
        const auto& helper = this->groupStateHelper();
        auto target = GuideRateModel::Target::OIL;
        for (std::string g = schedule.getWell(name, reportStepIdx).groupName(); !g.empty(); ) {
            const auto& grp = schedule.getGroup(g, reportStepIdx);
            if (grp.isProductionGroup()) {
                const auto cmode = grp.productionControls(summary_state).cmode;
                if (cmode != Group::ProductionCMode::NONE && cmode != Group::ProductionCMode::FLD) {
                    target = helper.getProductionGuideTargetModeFromControlMode(cmode);
                    break;
                }
            }
            if (grp.parent() == g) { break; }
            g = grp.parent();
        }
        return deckGuideOn(name, target);
    };

    std::map<std::string, Scalar> new_pressures;
    std::set<std::string> decided;
    std::set<std::string> route_wells;
    std::map<std::string, Well::ProducerCMode> assigned_cmode;
    std::map<std::string, std::vector<Scalar>> assigned_rates;
    std::map<std::string, bool> dead_now;
    // A well at zero rate given a rate again more often than the budget stays shut for the
    // rest of the step, as legacy's max_well_status_switch: a well on its lift cliff otherwise
    // flips every Newton iteration and the step never converges.
    const auto step_key = std::pair{static_cast<double>(simulator_.time()), static_cast<double>(simulator_.timeStepSize())};
    if (step_key != this->controller_revival_step_) {
        this->controller_revival_step_ = step_key;
        this->controller_revivals_.clear();
    }
    auto held_dead = [&](const std::string& name, const Scalar rate_now) {
        const auto r = this->controller_revivals_.find(name);
        return !(rate_now > Scalar{0}) && r != this->controller_revivals_.end()
            && r->second >= param_.group_controller_max_revivals_;
    };
    std::set<std::string> kept_dead;
    std::set<std::string> repaired_dead;   // shut by a repair of a rejected answer
    std::set<std::string> route_groups, switched;
    std::map<std::string, Group::ProductionCMode> holding;
    int n_judged_wrong = 0;

    struct RootSpec {
        std::string root_name;
        std::optional<Scalar> pressure;
        const std::string& name() const { return root_name; }
    };
    std::vector<RootSpec> roots;
    if (no_network) {
        roots.push_back({"FIELD", Scalar{0}});      // no node below it, its pressure is never used
    } else {
        for (const auto& root_ref : network.roots()) {
            const auto& r = root_ref.get();
            roots.push_back({r.name(), r.terminal_pressure().has_value()
                                           ? std::optional<Scalar>(static_cast<Scalar>(*r.terminal_pressure()))
                                           : std::nullopt});
        }
    }
    for (const auto& root : roots) {
        if (!root.pressure.has_value()) {
            return giveUp(fmt::format("{} has no terminal pressure", root.name()));
        }
        const Scalar terminal = *root.pressure;
        Sys system(*this->getVFPProperties().getProd(), units);
        system.setTerminalPressure(terminal);
        // Nodes, parents before children.
        std::map<std::string, int> index;
        std::vector<std::string> order{root.name()};
        system.addNode(NetworkSolve::Node{order.front(), -1, NetworkSolve::NoTable}, Scalar{0});
        index[order.front()] = 0;
        for (std::size_t at = 0; !no_network && at < order.size(); ++at) {
            for (const auto& branch : network.downtree_branches(order[at])) {
                const auto& child = branch.downtree_node();
                if (index.count(child)) {
                    continue;
                }
                index[child] = static_cast<int>(order.size());
                order.push_back(child);
                const auto& child_node = network.node(child);
                if (child_node.terminal_pressure().has_value()) {
                    return giveUp(fmt::format("{} is a fixed-pressure node below the root", child));
                }
                if (child_node.as_choke()) {
                    return giveUp(fmt::format("{} is an autochoke node", child));
                }
                Scalar alq = 0.0;
                if (branch.vfp_table().has_value()) {
                    const auto& table = this->getVFPProperties().getProd()->getTable(*branch.vfp_table());
                    alq = branch.alq_value(VFPProdTable::ALQDimension(table.getALQType(), units)).value_or(0.0);
                }
                NetworkSolve::Node node{child, static_cast<int>(at),
                                        branch.vfp_table().value_or(NetworkSolve::NoTable)};
                node.efficiency = child_node.efficiency();
                system.addNode(std::move(node), alq);
                // Satellite production: a rate with no well behind it.
                const auto& grp = schedule.getGroup(child, reportStepIdx);
                if (grp.hasSatelliteProduction()) {
                    const auto& gsat = schedule[reportStepIdx].satelliteProduction;
                    if (gsat.has(child)) {
                        const auto r = gsat.get(child).getRates(summary_state);
                        std::array<Scalar, Sys::NP> source{
                            static_cast<Scalar>(r[GSatProd::Rate::Water]),
                            static_cast<Scalar>(r[GSatProd::Rate::Oil]),
                            static_cast<Scalar>(r[GSatProd::Rate::Gas])};
                        if (child_node.add_gas_lift_gas()) {
                            source[2] += static_cast<Scalar>(r[GSatProd::Rate::GLift]);
                        }
                        system.setNodeSource(index.at(child), source);
                    }
                }
            }
        }
        // The wells under it.
        struct TreeWell { int index; std::string name; Scalar deck_guide; };
        std::vector<TreeWell> tree_wells;
        std::vector<std::string> own_control;     // in the system, free on their own limits
        std::map<std::string, Well::ProducerCMode> own_limit_mode;   // which own limit the allowance is
        // Every own rate limit as the oil rate at which it binds on the well's
        // linear IPR, and the smallest of them; exact for the line.
        auto allowance = [&](const typename Sys::Well& w, const Well::ProductionControls& c) {
            Scalar best = c.hasControl(Well::ProducerCMode::ORAT) && c.oil_rate > 0 ? static_cast<Scalar>(c.oil_rate) : Scalar{0};
            Well::ProducerCMode mode = Well::ProducerCMode::ORAT;
            auto consider = [&](const Well::ProducerCMode m, const Scalar limit, const int ph_a, const int ph_b) {
                if (!c.hasControl(m) || !(limit > 0)) return;
                // rate on the limit's phases q = A + B*bhp; bhp where it equals the limit; oil there.
                const Scalar A = w.ipr_a[ph_a] + (ph_b >= 0 ? w.ipr_a[ph_b] : Scalar{0});
                const Scalar B = w.ipr_b[ph_a] + (ph_b >= 0 ? w.ipr_b[ph_b] : Scalar{0});
                if (!(B < Scalar{0})) return;
                const Scalar bhp = (limit - A) / B;
                const Scalar oil = std::max(w.ipr_a[1] + w.ipr_b[1] * bhp, Scalar{0});
                if (!(best > Scalar{0}) || oil < best) { best = oil; mode = m; }
            };
            consider(Well::ProducerCMode::WRAT, static_cast<Scalar>(c.water_rate), 0, -1);
            consider(Well::ProducerCMode::GRAT, static_cast<Scalar>(c.gas_rate), 2, -1);
            consider(Well::ProducerCMode::LRAT, static_cast<Scalar>(c.liquid_rate), 0, 1);
            return std::pair{best, mode};
        };
        for (const auto& name : schedule.wellNames(reportStepIdx)) {
            const auto& well = schedule.getWell(name, reportStepIdx);
            if (!well.isProducer() || !well.predictionMode() || (!no_network && !index.count(well.groupName()))) {
                continue;
            }
            if (schedule.getGroup(well.groupName(), reportStepIdx).hasSatelliteProduction()) {
                continue;
            }
            if (schedule[reportStepIdx].glo().has_well(name)) {
                return giveUp(fmt::format("{} is under gas lift optimisation", name));
            }
            if (!this->wellState().has(name)) {
                continue;
            }
            const auto& ws = this->wellState().well(name);
            if (ws.status != WellStatus::OPEN) {
                continue;
            }
            const auto controls = well.productionControls(summary_state);
            const bool usable = static_cast<int>(ws.implicit_ipr_b.size()) >= pu.numActivePhases()
                && ws.implicit_ipr_b[pos[1]] > Scalar{0};
            const Scalar current = std::max(-ws.surface_rates[pos[1]], Scalar{0});
            const bool owned = groupControllable(well);
            typename Sys::Well w;
            w.name = name;
            w.group_controllable = owned;
            w.node = no_network ? 0 : index.at(well.groupName());
            w.vfp_table = controls.vfp_table_number;
            if (no_network && w.vfp_table > 0 && controls.thp_limit > 0.0) {
                w.own_thp = static_cast<Scalar>(controls.thp_limit);
            } else if (no_network) {
                w.vfp_table = 0;            // no thp limit: the table constrains nothing
            }
            if (w.vfp_table > 0) {
                // The table's datum is not the well's reference depth.
                const auto& wi = this->getWell(name);
                w.vfp_dp = wellhelpers::computeHydrostaticCorrection(
                    wi.refDepth(), this->getVFPProperties().getProd()->getTable(w.vfp_table).getDatumDepth(),
                    wi.refDensity(), wi.gravity());
                w.explicit_wfr = this->getVFPProperties().getExplicitWFR(w.vfp_table, wi.indexOfWell());
                w.explicit_gfr = this->getVFPProperties().getExplicitGFR(w.vfp_table, wi.indexOfWell());
                // The well model's flag, including the one its own failed solves set: ignoring
                // it costs 5-8 % oil on the network decks.
                w.explicit_vfp = wi.useVfpExplicit();
            }
            w.bhp_limit = static_cast<Scalar>(controls.bhp_limit);
            w.efficiency = static_cast<Scalar>(well.getEfficiencyFactor(/*network=*/true)) * ws.efficiency_scaling_factor;
            w.alq = ws.alq_state.get();
            w.node_adds_lift_gas = !no_network && network.node(well.groupName()).add_gas_lift_gas();
            if (w.node_adds_lift_gas) {
                w.lift_gas = w.alq;
            }
            const bool held = held_dead(name, current);
            if (held) {
                kept_dead.insert(name);
            }
            if (!usable || held) {
                if (!owned) {
                    continue;             // at zero rate on its own control; not part of it
                }
                // No inflow: shut in the system, counted so its group keeps its target.
                w.shut = true;
                w.guide = deckGuide(name);
                tree_wells.push_back({static_cast<int>(system.numWells()), name, w.guide});
                system.addWell(std::move(w));
                continue;
            }
            for (int ph = 0; ph < Sys::NP; ++ph) {
                // The well state holds q = b*bhp - a with production negative;
                // the system wants production positive, falling with bhp.
                w.ipr_a[ph] = ws.implicit_ipr_a[pos[ph]];
                w.ipr_b[ph] = -ws.implicit_ipr_b[pos[ph]];
            }
            if (current > Scalar{0}) {
                // Half the ratio the well has now: see Sys::ipr().
                for (const int ph : {0, 2}) {
                    w.min_ratio[ph] = Scalar{0.5} * std::max(-ws.surface_rates[pos[ph]], Scalar{0}) / current;
                }
            }
            w.q_start = current;
            if (owned) {
                const auto [allow, mode] = allowance(w, controls);
                w.oil_rate_limit = allow;
                own_limit_mode[name] = mode;
                w.guide = current;
                tree_wells.push_back({static_cast<int>(system.numWells()), name, deckGuide(name)});
            } else if (w.vfp_table > 0) {
                // On THP or on an own limit it may not reach at this node pressure: the
                // network places it. Pinned at its current rate it could not switch to THP.
                const auto [allow, mode] = allowance(w, controls);
                w.oil_rate_limit = allow;
                own_limit_mode[name] = mode;
                w.guide = std::max(current, w.oil_rate_limit);
                own_control.push_back(name);
            } else {
                if (!(current > Scalar{0})) {
                    continue;
                }
                w.oil_rate_limit = current;   // a source at what it does now
                w.guide = current;
                w.pinned = true;
            }
            system.addWell(std::move(w));
        }
        if (system.numWells() == 0) {
            return giveUp(fmt::format("no producers under {}", root.name()));
        }
        // Every deck limit of every group in the tree, and the satellite production counted
        // on it: the route holds one limit per group, chosen again from its own answer.
        std::map<std::string, int> gidx;
        std::map<int, std::vector<std::pair<Mode, Scalar>>> group_limits;
        std::map<int, std::array<Scalar, Sys::NP>> satellite_on;
        std::function<void()> refreshGuides = [] {};
        // The deck's group tree, parents first.
        if (!tree_wells.empty()) {
            const bool deck_guides = std::all_of(tree_wells.begin(), tree_wells.end(),
                                                 [](const TreeWell& t) { return t.deck_guide > Scalar{0}; });
            if (deck_guides) {
                for (const auto& t : tree_wells) {
                    system.setWellGuide(t.index, t.deck_guide);
                }
            }
            std::map<std::string, Scalar> wguide;
            std::set<std::string> here;
            for (const auto& t : tree_wells) {
                wguide[t.name] = system.wells()[t.index].guide;
            }
            // Every well in the system counts as here, including those not available for group
            // control: they are not held to a share, but their rates are in the group's total.
            for (const auto& w : system.wells()) { here.insert(w.name); }
            std::function<Scalar(const std::string&)> subtreeGuide = [&](const std::string& g) {
                const auto& grp = schedule.getGroup(g, reportStepIdx);
                Scalar sum = 0;
                for (const auto& child : grp.groups()) { sum += subtreeGuide(child); }
                for (const auto& wn : grp.wells()) {
                    if (const auto it = wguide.find(wn); it != wguide.end()) { sum += it->second; }
                }
                return sum;
            };
            std::function<bool(const std::string&)> allHere = [&](const std::string& g) {
                const auto& grp = schedule.getGroup(g, reportStepIdx);
                for (const auto& child : grp.groups()) { if (!allHere(child)) { return false; } }
                for (const auto& wn : grp.wells()) {
                    const auto& well = schedule.getWell(wn, reportStepIdx);
                    // A well the well model has stopped produces nothing; the target
                    // applies to the rest. Only a producing well outside the system
                    // (another network, another rank) takes the target away.
                    const bool open = this->wellState().has(wn)
                        && this->wellState().well(wn).status == WellStatus::OPEN;
                    if (well.isProducer() && well.predictionMode() && open && !here.count(wn)) {
                        return false;
                    }
                }
                return true;
            };
            std::vector<std::pair<Mode, Scalar>> pending_limits;
            std::function<void(const std::string&, int)> addTree = [&](const std::string& g, const int parent) {
                const auto& grp = schedule.getGroup(g, reportStepIdx);
                typename Sys::Group node;
                node.name = g;
                node.parent = parent;
                node.efficiency = grp.getGroupEfficiencyFactor(/*network=*/true);
                node.guide = subtreeGuide(g);
                node.available = grp.productionGroupControlAvailable();
                if (grp.isProductionGroup() && allHere(g)) {
                    const auto ctl = grp.productionControls(summary_state);
                    using C = Group::ProductionCMode;
                    // One target per group in the system: of the deck's limits, the one
                    // most binding at the group's current rates. The controller's group
                    // check re-decides when another one takes over.
                    const auto& gs = this->groupState();
                    std::array<Scalar, Sys::NP> now{};   // water, oil, gas
                    if (gs.has_production_rates(g)) {
                        const auto& r = gs.production_rates(g);
                        for (int ph = 0; ph < Sys::NP; ++ph) { now[ph] = pos[ph] < static_cast<int>(r.size()) ? r[pos[ph]] : Scalar{0}; }
                    }
                    const auto& act = ctl.group_limit_action;
                    const bool all_rate = act.allRates == Group::ExceedAction::RATE;
                    Scalar worst = -1;
                    std::vector<std::pair<Mode, Scalar>> limits_here;
                    auto consider = [&](const C m, const Group::ExceedAction a, const Scalar limit, const Mode mode, const Scalar rate) {
                        if (!grp.has_control(m) || !(all_rate || a == Group::ExceedAction::RATE)) return;
                        const Scalar ratio = limit > Scalar{0} ? rate / limit : (m == ctl.cmode ? Scalar{1e30} : Scalar{-1});
                        const Scalar score = m == ctl.cmode ? std::max(ratio, Scalar{0}) + Scalar{1e-9} : ratio;
                        if (score > worst) { worst = score; node.mode = mode; node.target = limit; }
                        if (limit > Scalar{0}) { limits_here.emplace_back(mode, limit); }
                    };
                    consider(C::ORAT, act.oil, ctl.oil_target, Mode::Oil, now[1]);
                    consider(C::WRAT, act.water, ctl.water_target, Mode::Water, now[0]);
                    consider(C::GRAT, act.gas, ctl.gas_target, Mode::Gas, now[2]);
                    consider(C::LRAT, act.liquid, ctl.liquid_target, Mode::Liquid, now[0] + now[1]);
                    // An explicit zero target means produce nothing; the system reads zero as none.
                    if (!(node.target > Scalar{0}) && grp.has_control(ctl.cmode)
                        && (ctl.cmode == C::ORAT || ctl.cmode == C::LRAT || ctl.cmode == C::GRAT || ctl.cmode == C::WRAT)) {
                        node.target = Scalar{1e-12};
                    }
                    pending_limits = std::move(limits_here);
                }
                const int me = system.addGroup(std::move(node));
                if (pending_limits.size() > 1) { group_limits[me] = std::move(pending_limits); }
                pending_limits.clear();
                gidx[g] = me;
                for (const auto& child : grp.groups()) { addTree(child, me); }
            };
            addTree("FIELD", -1);
            // Satellite production reduces every target above it.
            {
                const auto& gs = this->groupState();
                const auto& groups = system.groups();
                for (int g = 0; g < system.numGroups(); ++g) {
                    if (!schedule.getGroup(groups[g].name, reportStepIdx).hasSatelliteProduction()
                        || !gs.has_production_rates(groups[g].name)) {
                        continue;
                    }
                    const auto& r = gs.production_rates(groups[g].name);
                    std::array<Scalar, Sys::NP> q{};
                    for (int ph = 0; ph < Sys::NP; ++ph) {
                        q[ph] = (pos[ph] < static_cast<int>(r.size())) ? r[pos[ph]] : Scalar{0};
                    }
                    Scalar eff = groups[g].efficiency;
                    for (int a = groups[g].parent; a >= 0; a = groups[a].parent) {
                        for (int ph = 0; ph < Sys::NP; ++ph) { satellite_on[a][ph] += eff * q[ph]; }
                        if (groups[a].target > Scalar{0}) {
                            const auto c = Sys::modeWeights(groups[a].mode, groups[a].resv_coeff);
                            Scalar on = 0;
                            for (int ph = 0; ph < Sys::NP; ++ph) { on += c[ph] * q[ph]; }
                            system.setGroupLimit(a, groups[a].mode, std::max(groups[a].target - eff * on, Scalar{1e-12}));
                        }
                        eff *= groups[a].efficiency;
                    }
                }
            }
            // Every well counts in its group's rate; the tree decides which can be held (WGRUPCON).
            for (int w = 0; w < system.numWells(); ++w) {
                const auto it = gidx.find(schedule.getWell(system.wells()[w].name, reportStepIdx).groupName());
                if (it != gidx.end()) { system.setWellGroup(w, it->second); }
            }
            // A share is measured on the mode of the group that hands it out, so that is
            // the phase a node's guide rate is wanted on: a well's on its group's mode, a
            // group's on its parent's. Done again whenever a group changes the limit it holds.
            std::map<int, std::vector<std::string>> below;
            for (const auto& t : tree_wells) {
                for (int a = gidx.at(schedule.getWell(t.name, reportStepIdx).groupName()); a >= 0;
                     a = system.groups()[a].parent) {
                    below[a].push_back(t.name);
                }
            }
            refreshGuides = [&, below, deck_guides]() {
                if (!deck_guides) {
                    return;
                }
                const auto& groups = system.groups();
                for (const auto& t : tree_wells) {
                    const int g = gidx.at(schedule.getWell(t.name, reportStepIdx).groupName());
                    system.setWellGuide(t.index, deckGuideOn(t.name, targetOf(groups[g].mode)));
                }
                for (int g = 0; g < system.numGroups(); ++g) {
                    const auto it = below.find(g);
                    if (groups[g].parent < 0 || it == below.end()) {
                        continue;
                    }
                    Scalar sum = 0;
                    for (const auto& name : it->second) {
                        sum += deckGuideOn(name, targetOf(groups[groups[g].parent].mode));
                    }
                    system.setGroupGuide(g, sum);
                }
            };
            refreshGuides();
            system.setGroupTree(true);
            system.setGroupActiveSet(true);
        }
        system.setAnalyticJacobian(true);
        system.setComplementarity(false);
        system.finish();
        if (!tree_wells.empty()) {
            system.finishGroups();
        }
        // From the node pressures the network last had.
        std::vector<Scalar> guess(order.size(), terminal);
        const auto& previous = this->network_.nodePressures();
        for (std::size_t n = 0; n < order.size(); ++n) {
            if (const auto it = previous.find(order[n]); it != previous.end() && it->second > Scalar{0}) {
                guess[n] = it->second;
            }
        }
        // Who decides the group tree's active set inside the route. OPM_CONTROLLER_GROUP_SET:
        //   stein (default)  the balancer's tree, at every evaluation, in place of the route's walk
        //   greedy           the route's own walk, and a greedy switch between a group's limits
        //   stein-limits     the route's own walk, the balancer's tree choosing between the limits
        static const std::string group_set = [] {
            const char* v = std::getenv("OPM_CONTROLLER_GROUP_SET");
            return std::string(v != nullptr ? v : "stein");
        }();
        const bool stein_set = group_set == "stein" && !tree_wells.empty();
        int stein_mode_switches = 0;
        ProdGroupTreeBalancer::Tree<Scalar> stein_last;   // the tree behind the set in force
        if (stein_set) {
            system.setLiftTolerance(param_.group_controller_network_tolerance_);
            system.setTreeAllocator([&](const std::vector<Scalar>& oil_capacity,
                                        const std::vector<char>& not_holdable,
                                        typename Sys::TreeDecision& d) -> bool {
                const auto& wells = system.wells();
                std::set<std::string> individual;
                for (int w = 0; w < system.numWells(); ++w) {
                    if (not_holdable[w]) { individual.insert(wells[w].name); }
                }
                auto phases = [](const typename Sys::Well& well, const Scalar q_oil) {
                    const Scalar bhp = (q_oil - well.ipr_a[1]) / well.ipr_b[1];
                    std::array<Scalar, 3> q{};     // oil, water, gas
                    q[0] = std::max(Sys::ipr(well, 1, bhp), Scalar{0});
                    q[1] = std::max(Sys::ipr(well, 0, bhp), Scalar{0});
                    q[2] = std::max(Sys::ipr(well, 2, bhp), Scalar{0});
                    return q;
                };
                std::unordered_map<std::string, std::pair<int, Scalar>> capacity;
                std::map<std::string, std::array<Scalar, 3>> rates;
                for (int w = 0; w < system.numWells(); ++w) {
                    if (!(oil_capacity[w] > Scalar{0}) || !(wells[w].ipr_b[1] < Scalar{0})) {
                        continue;
                    }
                    const auto at_cap = phases(wells[w], oil_capacity[w]);
                    capacity[wells[w].name] = {static_cast<int>(Well::ProducerCMode::THP),
                                               at_cap[0] + at_cap[1] + at_cap[2]};
                    rates[wells[w].name] = at_cap;
                }
                if (capacity.empty()) {
                    return false;
                }
                // The tree works at fixed phase fractions and the wells sit on their inflow
                // lines: a few sweeps put the fractions where the allocation puts the wells.
                DeferredLogger quiet;
                ProdGroupTreeBalancer::Tree<Scalar> tree;
                for (int sweep = 0; sweep < 6; ++sweep) {
                    bool valid = false;
                    tree = ProdGroupTreeBalancer::decideTree(
                        static_cast<const BlackoilWellModelGeneric<Scalar, IndexTraits>&>(*this), summary_state,
                        reportStepIdx, static_cast<Scalar>(param_.group_tree_balancer_tolerance_), capacity,
                        rates, quiet, valid, &individual);
                    if (!valid) {
                        // Say why, once more with a logger that is kept.
                        std::string what;
                        for (const auto& [name, r] : rates) {
                            what += fmt::format(" {}: capacity sum {:.4g}, at o/w/g {:.4g}/{:.4g}/{:.4g};", name,
                                                capacity.at(name).second * 86400.0, r[0] * 86400.0, r[1] * 86400.0,
                                                r[2] * 86400.0);
                        }
                        deferred_logger.debug(fmt::format("Controller: the balancer's tree is not valid (sweep {}):{}",
                                                          sweep, what));
                        bool again = false;
                        const auto bad = ProdGroupTreeBalancer::decideTree(
                            static_cast<const BlackoilWellModelGeneric<Scalar, IndexTraits>&>(*this), summary_state,
                            reportStepIdx, static_cast<Scalar>(param_.group_tree_balancer_tolerance_), capacity,
                            rates, deferred_logger, again, &individual);
                        for (const auto& [name, node] : bad) {
                            deferred_logger.debug(fmt::format(
                                "Controller:   {} {} category {} mode {} o/w/g {:.4g}/{:.4g}/{:.4g} limits{}", name,
                                node.type == ProdNodeType::Well ? "well" : "group", static_cast<int>(node.modeCategory),
                                static_cast<int>(node.mode), -node.rates[0] * 86400.0, -node.rates[1] * 86400.0,
                                -node.rates[2] * 86400.0, [&node] {
                                    std::string l;
                                    for (const auto& [m, v] : node.Limits) {
                                        l += fmt::format(" {}={:.4g}", static_cast<int>(m), v * 86400.0);
                                    }
                                    return l;
                                }()));
                        }
                        return false;
                    }
                    Scalar moved = 0;
                    for (int w = 0; w < system.numWells(); ++w) {
                        const auto it = tree.find(wells[w].name);
                        if (it == tree.end() || capacity.count(wells[w].name) == 0) {
                            continue;
                        }
                        const Scalar oil = std::clamp(-it->second.rates[0], Scalar{0}, oil_capacity[w]);
                        const auto at = phases(wells[w], oil > Scalar{0} ? oil : oil_capacity[w]);
                        auto& r = rates[wells[w].name];
                        for (int c = 0; c < 3; ++c) {
                            moved = std::max(moved, std::abs(at[c] - r[c]) / std::max(r[c], Scalar{1e-12}));
                        }
                        r = at;
                    }
                    if (moved < Scalar{1e-3}) {
                        break;
                    }
                }
                // A group on a limit the system does not hold for it: hold that one.
                bool switched = false;
                for (const auto& [g, limits] : group_limits) {
                    const auto& grp = system.groups()[g];
                    const auto it = tree.find(grp.name);
                    if (it == tree.end() || it->second.modeCategory != ProdNodeModeCategory::Individual) {
                        continue;
                    }
                    Mode held = Mode::None;
                    switch (it->second.mode) {
                    case Well::ProducerCMode::ORAT: held = Mode::Oil; break;
                    case Well::ProducerCMode::WRAT: held = Mode::Water; break;
                    case Well::ProducerCMode::GRAT: held = Mode::Gas; break;
                    case Well::ProducerCMode::LRAT: held = Mode::Liquid; break;
                    default: break;
                    }
                    if (held == Mode::None || held == grp.mode) {
                        continue;
                    }
                    for (const auto& lim : limits) {
                        if (lim.first != held) {
                            continue;
                        }
                        const auto c = Sys::modeWeights(held, grp.resv_coeff);
                        Scalar sat = 0;
                        if (const auto is = satellite_on.find(g); is != satellite_on.end()) {
                            for (int ph = 0; ph < Sys::NP; ++ph) { sat += c[ph] * is->second[ph]; }
                        }
                        system.setGroupLimit(g, held, std::max(lim.second - sat, Scalar{1e-12}));
                        switched = true;
                        ++stein_mode_switches;
                    }
                }
                if (switched) {
                    refreshGuides();
                }
                d.bind.assign(system.numGroups(), Sys::GroupBind::Free);
                d.held.assign(system.numWells(), 0);
                d.oil.assign(system.numWells(), Scalar{0});
                for (int w = 0; w < system.numWells(); ++w) {
                    const auto it = tree.find(wells[w].name);
                    if (it == tree.end() || wells[w].group < 0) {
                        continue;
                    }
                    if (it->second.modeCategory == ProdNodeModeCategory::Group) {
                        d.held[w] = 1;
                        d.oil[w] = std::max(-it->second.rates[0], Scalar{0});
                    }
                }
                // The groups the tree holds on a limit of their own; every group between a
                // held well and the nearest such group above it passes that share on.
                std::vector<char> own(system.numGroups(), 0), carries(system.numGroups(), 0);
                for (int g = 0; g < system.numGroups(); ++g) {
                    const auto it = tree.find(system.groups()[g].name);
                    own[g] = it != tree.end() && system.groups()[g].target > Scalar{0}
                        && it->second.modeCategory == ProdNodeModeCategory::Individual;
                }
                for (int w = 0; w < system.numWells(); ++w) {
                    if (!d.held[w]) {
                        continue;
                    }
                    int top = wells[w].group;
                    while (top >= 0 && !own[top]) {
                        top = system.groups()[top].parent;
                    }
                    if (top < 0) {
                        d.held[w] = 0;        // nothing above it holds a limit
                        continue;
                    }
                    for (int g = wells[w].group; g != top; g = system.groups()[g].parent) {
                        carries[g] = 1;
                    }
                    carries[top] = 1;
                }
                for (int g = 0; g < system.numGroups(); ++g) {
                    if (carries[g]) {
                        d.bind[g] = own[g] ? Sys::GroupBind::Own : Sys::GroupBind::Share;
                    }
                }
                if (std::getenv("OPM_CONTROLLER_TRACE") != nullptr) {
                    std::string t = fmt::format("CTRLTRACE step={} it={} allocator individual {{", reportStepIdx,
                                                simulator_.problem().iterationContext().iteration());
                    for (const auto& n : individual) { t += " " + n; }
                    t += " }";
                    for (const auto& [name, node] : tree) {
                        t += fmt::format(" | {} {} cat {} mode {} oil {:.1f} wat {:.1f}", name,
                                         node.type == ProdNodeType::Well ? "w" : "g",
                                         static_cast<int>(node.modeCategory), static_cast<int>(node.mode),
                                         -node.rates[0] * 86400.0, -node.rates[1] * 86400.0);
                    }
                    deferred_logger.debug(t);
                }
                stein_last = std::move(tree);
                return true;
            });
        }
        // Diagnostic: OPM_CONTROLLER_CLIFF=hold keeps a well at the rate it had before a cliff
        // instead of letting it die inside the route's solve.
        static const auto cliff_rule = [] {
            const char* v = std::getenv("OPM_CONTROLLER_CLIFF");
            return (v != nullptr && std::string(v) == "hold") ? NetworkSolve::CliffRule::Hold
                                                               : NetworkSolve::CliffRule::Die;
        }();
        const NetworkSolve::Parameters<Scalar> params{Scalar{1e-2}, 50};
        // Stein's continuation: every inflow line crosses its tubing curve once, the shut decision
        // is taken after the solve for every well at once.
        const bool extension = param_.group_controller_tubing_extension_;
        auto solveRoute = [&](const std::vector<Scalar>& start, const bool keep_dead) {
            if (!extension) {
                return NetworkSolve::solveReduced(system, start, params, /*eliminate=*/true, cliff_rule, keep_dead);
            }
            auto ex = NetworkSolve::solveReducedOnExtension(system, start, params, NetworkSolve::Closing::All,
                                                            20, keep_dead);
            auto r = std::move(ex.last);
            r.converged = r.converged && ex.converged;
            r.iterations = ex.iterations;
            r.evaluations = ex.evaluations;
            return r;
        };
        auto rr = solveRoute(guess, false);
        // A shut decided on the inflow line alone is confirmed with the well model: a line taken far
        // from the operating point (a well throttled by its group) can miss the tubing where the well
        // lifts fine. Where the well model finds a crossing at the well-head pressure, the lines are
        // re-anchored through it and the route solves again.
        // Once per well: at a network cliff the well lifts at the pressure of the shut answer and still
        // has no fixed point flowing; if the re-solve shuts it again, that shut stands.
        std::set<int> confirmed;
        for (int round = 0; rr.converged && round < param_.group_controller_max_repairs_; ++round) {
            std::string revived;
            for (int w = 0; w < system.numWells(); ++w) {
                const auto& well = system.wells()[w];
                if (system.controlLetter(w) != 'S' || well.shut || well.pinned || well.vfp_table <= 0
                    || !confirmed.insert(w).second) { continue; }
                const auto wit = std::find_if(well_container_.begin(), well_container_.end(),
                                              [&](const auto& x) { return x->name() == well.name; });
                if (wit == well_container_.end()) { continue; }
                auto& wi = *wit;
                const Scalar p_w = well.own_thp > Scalar{0} ? well.own_thp
                                 : well.node == 0 ? terminal : rr.node_pressure[well.node];
                const auto saved = wi->getDynamicThpLimit();
                wi->setDynamicThpLimit(p_w);
                const auto bhp = wi->computeBhpAtThpLimitProdWithAlq(simulator_, this->groupStateHelper(),
                                                                     summary_state, well.alq, false);
                wi->setDynamicThpLimit(saved);
                if (!bhp) { continue; }
                std::vector<Scalar> q(pu.numActivePhases(), Scalar{0});
                wi->computeWellRatesWithBhp(simulator_, *bhp, q, deferred_logger);
                std::array<Scalar, Sys::NP> a = well.ipr_a;
                bool flows = false;
                for (int ph = 0; ph < Sys::NP; ++ph) {
                    const Scalar q_true = pos[ph] >= 0 ? std::max(-q[pos[ph]], Scalar{0}) : Scalar{0};
                    a[ph] = q_true - well.ipr_b[ph] * *bhp;
                    flows = flows || (ph == 1 && q_true > Scalar{0});
                }
                if (!flows) { continue; }
                system.setWellIpr(w, a, well.ipr_b);
                system.setWellDeadAbove(w, Scalar{0});
                system.reviveWell(w, Sys::ipr(system.wells()[w], 1, *bhp));
                revived += fmt::format(" {} ({:.1f} sm3/d at {:.2f} bar)", well.name,
                                       Sys::ipr(system.wells()[w], 1, *bhp) * 86400.0, p_w / 1e5);
            }
            if (revived.empty()) { break; }
            deferred_logger.debug(fmt::format("Controller: the well model lifts wells the route shut under {} at report step "
                                              "{}:{}; solved again", root.name(), reportStepIdx, revived));
            rr = solveRoute(rr.node_pressure, true);
        }
        // A group with several limits holds the one its own answer violates most: solve,
        // look at every limit, switch and solve again. Bounded; a repeat ends it.
        // OPM_CONTROLLER_STEIN_LIMITS: the balancer's tree on the route's answer picks the
        // limit each group holds, in place of the greedy switch below.
        const bool stein_limits = group_set == "stein-limits";
        int limit_switches = 0;
        Scalar stein_gap = 0;
        std::string stein_gap_well;
        bool stein_valid = true;
        for (int round = 0; !stein_set && rr.converged && round < 4 && !group_limits.empty(); ++round) {
            if (stein_limits) {
                // Every well's capacity and phase split where the route put it.
                std::unordered_map<std::string, std::pair<int, Scalar>> capacity;
                std::map<std::string, std::array<Scalar, 3>> at_answer;
                for (int w = 0; w < system.numWells(); ++w) {
                    const auto& well = system.wells()[w];
                    if (!(well.ipr_b[1] < Scalar{0})) {
                        continue;
                    }
                    auto phases = [&well](const Scalar q_oil) {
                        const Scalar bhp = (q_oil - well.ipr_a[1]) / well.ipr_b[1];
                        std::array<Scalar, 3> q{};     // oil, water, gas
                        q[0] = std::max(Sys::ipr(well, 1, bhp), Scalar{0});
                        q[1] = std::max(Sys::ipr(well, 0, bhp), Scalar{0});
                        q[2] = std::max(Sys::ipr(well, 2, bhp), Scalar{0});
                        return q;
                    };
                    Scalar cap = well.ipr_a[1] + well.ipr_b[1] * well.bhp_limit;
                    auto mode = Well::ProducerCMode::BHP;
                    if (well.pinned) {
                        cap = well.oil_rate_limit;
                    } else if (well.vfp_table > 0) {
                        const Scalar p_node = well.own_thp > Scalar{0} ? well.own_thp
                            : well.node == 0 ? terminal : rr.node_pressure[well.node];
                        const Scalar lift = system.thpPotential(well, p_node);
                        if (lift < cap) {
                            cap = lift;
                            mode = Well::ProducerCMode::THP;
                        }
                    }
                    if (!(cap > Scalar{0})) {
                        continue;
                    }
                    const auto at_cap = phases(cap);
                    capacity[well.name] = {static_cast<int>(mode), at_cap[0] + at_cap[1] + at_cap[2]};
                    at_answer[well.name] = rr.well_rate[w] > Scalar{0} ? phases(std::min(rr.well_rate[w], cap)) : at_cap;
                }
                bool valid = false;
                const auto tree = ProdGroupTreeBalancer::decideTree(
                    static_cast<const BlackoilWellModelGeneric<Scalar, IndexTraits>&>(*this), summary_state,
                    reportStepIdx, static_cast<Scalar>(param_.group_tree_balancer_tolerance_), capacity,
                    at_answer, deferred_logger, valid);
                stein_valid = valid;
                stein_gap = 0;
                for (int w = 0; w < system.numWells(); ++w) {
                    const auto it = tree.find(system.wells()[w].name);
                    if (it == tree.end()) {
                        continue;
                    }
                    const Scalar mine = rr.well_rate[w], theirs = -it->second.rates[0];
                    const Scalar gap = std::abs(mine - theirs) / std::max({mine, theirs, Scalar{1e-9}});
                    if (gap > stein_gap) {
                        stein_gap = gap;
                        stein_gap_well = system.wells()[w].name;
                    }
                }
                bool switched = false;
                for (const auto& [g, limits] : group_limits) {
                    const auto& grp = system.groups()[g];
                    const auto it = tree.find(grp.name);
                    if (it == tree.end() || it->second.modeCategory != ProdNodeModeCategory::Individual) {
                        continue;
                    }
                    Mode held = Mode::None;
                    switch (it->second.mode) {
                    case Well::ProducerCMode::ORAT: held = Mode::Oil; break;
                    case Well::ProducerCMode::WRAT: held = Mode::Water; break;
                    case Well::ProducerCMode::GRAT: held = Mode::Gas; break;
                    case Well::ProducerCMode::LRAT: held = Mode::Liquid; break;
                    default: break;
                    }
                    if (held == Mode::None || held == grp.mode) {
                        continue;
                    }
                    for (const auto& lim : limits) {
                        if (lim.first != held) {
                            continue;
                        }
                        const auto c = Sys::modeWeights(held, grp.resv_coeff);
                        Scalar sat = 0;
                        if (const auto is = satellite_on.find(g); is != satellite_on.end()) {
                            for (int ph = 0; ph < Sys::NP; ++ph) { sat += c[ph] * is->second[ph]; }
                        }
                        system.setGroupLimit(g, held, std::max(lim.second - sat, Scalar{1e-12}));
                        switched = true;
                        ++limit_switches;
                    }
                }
                if (!switched) {
                    break;
                }
                refreshGuides();
                rr = NetworkSolve::solveReduced(system, rr.node_pressure, params, /*eliminate=*/true,
                                                cliff_rule);
                continue;
            }
            std::map<int, std::array<Scalar, Sys::NP>> produced = satellite_on;
            for (int w = 0; w < system.numWells(); ++w) {
                const auto& well = system.wells()[w];
                const Scalar q_oil = rr.well_rate[w];
                const auto gi = gidx.find(schedule.getWell(well.name, reportStepIdx).groupName());
                if (!(q_oil > Scalar{0}) || !(well.ipr_b[1] < Scalar{0}) || gi == gidx.end()) {
                    continue;
                }
                const Scalar bhp = (q_oil - well.ipr_a[1]) / well.ipr_b[1];
                Scalar eff = well.efficiency;
                for (int a = gi->second; a >= 0; a = system.groups()[a].parent) {
                    for (int ph = 0; ph < Sys::NP; ++ph) {
                        produced[a][ph] += eff * std::max(Sys::ipr(well, ph, bhp), Scalar{0});
                    }
                    eff *= system.groups()[a].efficiency;
                }
            }
            bool switched = false;
            for (const auto& [g, limits] : group_limits) {
                const auto& grp = system.groups()[g];
                Scalar worst = Scalar{1} + Scalar{1e-3};
                const std::pair<Mode, Scalar>* take = nullptr;
                for (const auto& lim : limits) {
                    const auto c = Sys::modeWeights(lim.first, grp.resv_coeff);
                    Scalar on = 0;
                    for (int ph = 0; ph < Sys::NP; ++ph) { on += c[ph] * produced[g][ph]; }
                    if (lim.first != grp.mode && on / lim.second > worst) {
                        worst = on / lim.second;
                        take = &lim;
                    }
                }
                if (take != nullptr) {
                    const auto c = Sys::modeWeights(take->first, grp.resv_coeff);
                    Scalar sat = 0;
                    if (const auto it = satellite_on.find(g); it != satellite_on.end()) {
                        for (int ph = 0; ph < Sys::NP; ++ph) { sat += c[ph] * it->second[ph]; }
                    }
                    system.setGroupLimit(g, take->first, std::max(take->second - sat, Scalar{1e-12}));
                    switched = true;
                    ++limit_switches;
                }
            }
            if (!switched) {
                break;
            }
            refreshGuides();
            rr = NetworkSolve::solveReduced(system, rr.node_pressure, params, /*eliminate=*/true,
                                            cliff_rule);
        }
        if (limit_switches > 0) {
            deferred_logger.debug(fmt::format("Controller: {} group limit switches inside the decision under {}",
                                              limit_switches, root.name()));
        }
        if (stein_set) {
            deferred_logger.debug(fmt::format("Controller: the balancer's tree set the groups under {}: {} calls, "
                                              "{} fell back to the route's walk, {} limit switches", root.name(),
                                              system.treeAllocatorCalls(), system.treeAllocatorFallbacks(),
                                              stein_mode_switches));
        }
        if (stein_limits && !group_limits.empty()) {
            deferred_logger.debug(fmt::format("Controller: balancer on the route's answer under {}: {}, "
                                              "allocation differs by at most {:.1f} % ({})", root.name(),
                                              stein_valid ? "valid" : "NOT valid", 100.0 * stein_gap, stein_gap_well));
        }
        if (!rr.converged) {
            return giveUp(fmt::format("the reduced route did not converge under {} in {} iterations, residual {:.3g}",
                                      root.name(), rr.iterations, rr.residual));
        }
        // The judge: every row of the route's own problem, before anything is written.
        std::string letters;
        auto judge = [&] {
            letters.clear();
            for (int w = 0; w < system.numWells(); ++w) { letters += system.controlLetter(w); }
            return NetworkSolve::verifyAnswer(system, rr.node_pressure, rr.well_rate, letters,
                                              param_.group_controller_network_tolerance_,
                                              param_.group_controller_rate_tolerance_);
        };
        auto verdict = judge();
        auto describe = [](const NetworkSolve::Verdict& v) {
            std::string what;
            for (const auto& [k, n] : v.violations) { what += fmt::format(" {}x{}", n, k); }
            return what;
        };
        // A rejected answer is repaired toward a solution of the route's problem with shut allowed:
        // a well that has no operating point is shut, and the problem solved again.
        for (int repair = 0; !verdict.ok && repair < param_.group_controller_max_repairs_; ++repair) {
            ++n_judged_wrong;
            std::vector<int> shut;
            std::string why;
            if (rr.on_cliff && verdict.violations.count("node pressure off its branch")) {
                // No fixed point with them flowing; dead is the consistent answer.
                shut = rr.cliff_wells;
                why = "no fixed point flowing";
            }
            // Held at a rate its tubing cannot lift: not holdable, and allocate again.
            std::vector<int> unhold;
            if (shut.empty()) {
                for (const int w : verdict.wells_unliftable) {
                    if (system.controlLetter(w) == 'R' || system.controlLetter(w) == 'G') { unhold.push_back(w); }
                }
            }
            std::string names;
            for (const int w : shut) { names += " " + system.wells()[w].name; }
            for (const int w : unhold) { names += " " + system.wells()[w].name; }
            if (shut.empty()) {
                for (const int w : verdict.wells_unliftable) {
                    names += fmt::format(" [cannot lift: {} {} {:.1f}]", system.wells()[w].name, system.controlLetter(w),
                                         rr.well_rate[w] * 86400.0);
                }
                const auto& groups = system.groups();
                for (std::size_t k = 0; k < verdict.groups_over.size(); ++k) {
                    const int g = verdict.groups_over[k];
                    names += fmt::format(" [{} mode {} target {:.1f} at {:.4f} of it:", groups[g].name,
                                         static_cast<int>(groups[g].mode), groups[g].target * 86400.0,
                                         verdict.over_ratio[k]);
                    for (int w = 0; w < system.numWells(); ++w) {
                        for (int a = system.wells()[w].group; a >= 0; a = groups[a].parent) {
                            if (a == g) {
                                names += fmt::format(" {} {} {:.1f}{}{}", system.wells()[w].name, system.controlLetter(w),
                                                     rr.well_rate[w] * 86400.0,
                                                     system.wells()[w].group_controllable ? "" : " (not group controllable)",
                                                     system.holdable(w) ? " holdable" : "");
                                break;
                            }
                        }
                    }
                    names += "]";
                }
            }
            deferred_logger.debug(fmt::format("Controller: the judge rejects the network answer under {} at report step {}:{}; "
                                              "{}", root.name(), reportStepIdx, describe(verdict),
                                              !shut.empty() ? "shut" + names + " (" + why + ")"
                                              : !unhold.empty() ? "not holdable" + names : "no repair" + names));
            if (shut.empty() && unhold.empty()) { break; }
            for (const int w : unhold) { system.forceNotHoldable(w); }
            for (const int w : shut) {
                system.shutWell(w);
                repaired_dead.insert(system.wells()[w].name);
            }
            rr = solveRoute(rr.node_pressure, true);
            if (!rr.converged) { break; }
            verdict = judge();
        }
        if (!verdict.ok || !rr.converged) {
            // Nothing written for this root: its wells keep the previous decision, and the hand-over
            // is marked unsettled.
            this->controller_rejected_ = true;
            deferred_logger.debug(fmt::format("Controller: under {} at report step {} no answer passes the judge ({}); "
                                              "the previous decision stands", root.name(), reportStepIdx,
                                              rr.converged ? describe(verdict) : std::string(" not converged")));
            for (const auto& w : system.wells()) {
                route_wells.insert(w.name);
                repaired_dead.erase(w.name);
                if (this->controller_decided_wells_.count(w.name)) {
                    decided.insert(w.name);
                    if (const auto it = this->controller_assigned_cmode_.find(w.name); it != this->controller_assigned_cmode_.end()) {
                        assigned_cmode[w.name] = it->second;
                    }
                    if (const auto it = this->controller_assigned_rates_.find(w.name); it != this->controller_assigned_rates_.end()) {
                        assigned_rates[w.name] = it->second;
                    }
                }
            }
            continue;
        }
        deferred_logger.debug(fmt::format("Controller: network route under {} at report step {}: {} iterations, {} evaluations, "
                                          "{} set changes, {} lookups, set {}, judge {}{}",
                                          root.name(), reportStepIdx, rr.iterations, rr.evaluations, rr.set_changes,
                                          system.lookups(), letters, verdict.ok ? "ok" : "REJECTS",
                                          rr.stalls || rr.cliffs ? fmt::format(" ({} stalls, {} cliffs)", rr.stalls, rr.cliffs) : ""));
        auto& st = this->controller_stats_;
        ++st.decisions;
        st.route_iterations += rr.iterations;
        st.route_evaluations += rr.evaluations;
        st.set_changes += rr.set_changes;
        st.lookups += system.lookups();
        for (std::size_t n = 0; n < order.size(); ++n) {
            new_pressures[order[n]] = rr.node_pressure[n];
        }
        // The route never touches the root's entry: it is the terminal pressure.
        new_pressures[order.front()] = terminal;
        if (no_network) {
            new_pressures.clear();      // there are no nodes; the wells keep their own thp limits
        }
        // OPM_CONTROLLER_TRACE: one line per well per decision; OPM_CONTROLLER_DUMP=<prefix>:
        // every system written for the bench to replay.
        static const bool trace = std::getenv("OPM_CONTROLLER_TRACE") != nullptr;
        static const char* dump_prefix = std::getenv("OPM_CONTROLLER_DUMP");
        if (trace) {
            deferred_logger.debug(fmt::format("CTRLTRACE step={} it={} tree={} groups={} wells={} letters={}",
                reportStepIdx, simulator_.problem().iterationContext().iteration(),
                system.usesGroupTree() ? 1 : 0, system.numGroups(), system.numWells(), letters));
            for (int g = 0; g < system.numGroups(); ++g) {
                const auto& grp = system.groups()[g];
                deferred_logger.debug(fmt::format("CTRLTRACE step={} it={} group {} parent={} target={:.4g} mode={} guide={:.4g} avail={}",
                    reportStepIdx, simulator_.problem().iterationContext().iteration(), grp.name, grp.parent,
                    grp.target * 86400.0, static_cast<int>(grp.mode), grp.guide * 86400.0, grp.available ? 1 : 0));
            }
            for (int w = 0; w < system.numWells(); ++w) {
                const auto& well = system.wells()[w];
                const Scalar p = well.own_thp > Scalar{0} ? well.own_thp : rr.node_pressure[well.node];
                deferred_logger.debug(fmt::format(
                    "CTRLTRACE step={} it={} {} {} q={:.1f} qstart={:.1f} limit={:.1f} p={:.3f} "
                    "ipr_oil=({:.4g},{:.4g}) thp_pot={:.1f} bhp_lim={:.2f}{} ipr_from='{}'",
                    reportStepIdx, simulator_.problem().iterationContext().iteration(), well.name,
                    system.controlLetter(w), rr.well_rate[w] * 86400.0, well.q_start * 86400.0,
                    well.oil_rate_limit * 86400.0, p / 1e5, well.ipr_a[1] * 86400.0, well.ipr_b[1] * 86400.0 * 1e5,
                    well.vfp_table > 0 ? system.thpPotential(well, p) * 86400.0 : -1.0, well.bhp_limit / 1e5,
                    well.pinned ? " pinned" : "",
                    this->controller_ipr_source_.count(well.name) ? this->controller_ipr_source_.at(well.name) : std::string("-")));
                // The gas side: the line, the system's gas at its oil answer, the well's current gas.
                const Scalar qo = rr.well_rate[w];
                const Scalar bhp_w = (qo > Scalar{0} && well.ipr_b[1] < Scalar{0}) ? (qo - well.ipr_a[1]) / well.ipr_b[1] : well.bhp_limit;
                const Scalar gas_sys = std::max(well.ipr_a[2] + well.ipr_b[2] * bhp_w, Scalar{0});
                const Scalar gas_now = this->wellState().has(well.name)
                    ? -this->wellState().well(well.name).surface_rates[pos[2]] : Scalar{0};
                if (well.vfp_table > 0 && this->wellState().has(well.name)) {
                    // The tubing as the route sees it, at the rates and bhp the well has now.
                    const auto& wsn = this->wellState().well(well.name);
                    std::array<Scalar, Sys::NP> now{};
                    for (int ph = 0; ph < Sys::NP; ++ph) { now[ph] = std::max(-wsn.surface_rates[pos[ph]], Scalar{0}); }
                    deferred_logger.debug(fmt::format(
                        "CTRLTRACE step={} it={} {} tubing: table {} alq {:.4g} vfp_dp {:.3f} bar{}; at the well's rates "
                        "w/o/g {:.1f}/{:.1f}/{:.0f} and thp {:.2f}: route bhp {:.2f}, well bhp {:.2f}, well thp {:.2f}",
                        reportStepIdx, simulator_.problem().iterationContext().iteration(), well.name, well.vfp_table,
                        well.alq, well.vfp_dp / 1e5,
                        well.explicit_vfp ? fmt::format(", explicit wfr {:.3g} gfr {:.4g}", well.explicit_wfr, well.explicit_gfr) : std::string(),
                        now[0] * 86400.0, now[1] * 86400.0, now[2] * 86400.0, p / 1e5,
                        (system.tubingBhp(well, p, now) - well.vfp_dp) / 1e5, wsn.bhp / 1e5, wsn.thp / 1e5));
                    // The well model's own answer at this thp, with its true inflow, and its true
                    // inflow at the route's bhp: where the route's line and the well part ways.
                    const auto wit = std::find_if(well_container_.begin(), well_container_.end(),
                                                  [&](const auto& x) { return x->name() == well.name; });
                    if (wit == well_container_.end()) { continue; }
                    auto& wi = *wit;
                    const auto saved = wi->getDynamicThpLimit();
                    wi->setDynamicThpLimit(p);
                    const auto own = wi->computeBhpAtThpLimitProdWithAlq(simulator_, this->groupStateHelper(),
                                                                         summary_state, well.alq, false);
                    // And at the node pressure the previous decision handed over, for a jump.
                    const auto& prev_p = this->network_.nodePressures();
                    const auto pp = prev_p.find(schedule.getWell(well.name, reportStepIdx).groupName());
                    std::optional<Scalar> own_prev;
                    const Scalar p_prev = pp != prev_p.end() ? pp->second : Scalar{0};
                    if (p_prev > Scalar{0} && std::abs(p_prev - p) > Scalar{0.1e5}) {
                        wi->setDynamicThpLimit(p_prev);
                        own_prev = wi->computeBhpAtThpLimitProdWithAlq(simulator_, this->groupStateHelper(),
                                                                       summary_state, well.alq, false);
                    }
                    wi->setDynamicThpLimit(saved);
                    if (p_prev > Scalar{0} && std::abs(p_prev - p) > Scalar{0.1e5}) {
                        std::vector<Scalar> pq(pu.numActivePhases(), Scalar{0});
                        if (own_prev) { wi->computeWellRatesWithBhp(simulator_, *own_prev, pq, deferred_logger); }
                        deferred_logger.debug(fmt::format(
                            "CTRLTRACE step={} it={} {} well model at the previous node pressure {:.2f}: {}",
                            reportStepIdx, simulator_.problem().iterationContext().iteration(), well.name,
                            p_prev / 1e5, own_prev ? fmt::format("bhp {:.2f}, oil {:.1f}", *own_prev / 1e5,
                                                                 -pq[pos[1]] * 86400.0)
                                                   : std::string("no crossing")));
                    }
                    std::vector<Scalar> tq(pu.numActivePhases(), Scalar{0});
                    wi->computeWellRatesWithBhp(simulator_, bhp_w, tq, deferred_logger);
                    std::vector<Scalar> oq(pu.numActivePhases(), Scalar{0});
                    if (own) { wi->computeWellRatesWithBhp(simulator_, *own, oq, deferred_logger); }
                    deferred_logger.debug(fmt::format(
                        "CTRLTRACE step={} it={} {} well model at thp {:.2f}: {}; true inflow at the route's bhp {:.2f}: "
                        "oil {:.1f} gas {:.0f} (route {:.1f} / {:.0f})",
                        reportStepIdx, simulator_.problem().iterationContext().iteration(), well.name, p / 1e5,
                        own ? fmt::format("bhp {:.2f}, oil {:.1f} gas {:.0f}", *own / 1e5, oq[pos[1]] * 86400.0,
                                          oq[pos[2]] * 86400.0)
                            : std::string("no crossing"),
                        bhp_w / 1e5, tq[pos[1]] * 86400.0, tq[pos[2]] * 86400.0, qo * 86400.0, gas_sys * 86400.0));
                }
                deferred_logger.debug(fmt::format("CTRLTRACE step={} it={} {} gas: ipr=({:.4g},{:.4g}) bhp={:.2f} gas_sys={:.0f} gas_now={:.0f} shut={}",
                    reportStepIdx, simulator_.problem().iterationContext().iteration(), well.name,
                    well.ipr_a[2] * 86400.0, well.ipr_b[2] * 86400.0 * 1e5, bhp_w / 1e5, gas_sys * 86400.0, gas_now * 86400.0, well.shut ? 1 : 0));
            }
        }
        if (dump_prefix) {
            std::ofstream out(fmt::format("{}_ctrl_{}.txt", dump_prefix, this->controller_dumps_written_++));
            if (out) { NetworkSolve::write(system, guess, out); }
        }
        // The controls and targets.
        const auto held = system.heldTargets(rr.well_rate);
        for (const auto& g : system.groups()) {
            route_groups.insert(g.name);
        }
        for (int w = 0; w < system.numWells(); ++w) {
            const auto& well = system.wells()[w];
            route_wells.insert(well.name);
            if (well.pinned || !this->wellState().has(well.name)) {
                continue;
            }
            auto& ws = this->wellState().well(well.name);
            const auto& t = held[w];
            using GC = Group::ProductionCMode;
            const auto cmode_before = ws.production_cmode;
            switch (t.control) {
            case Ctrl::Tree:
            case Ctrl::Grup: {
                if (t.group < 0) {
                    break;
                }
                const GC cmode = t.mode == Mode::Gas    ? GC::GRAT
                               : t.mode == Mode::Water  ? GC::WRAT
                               : t.mode == Mode::Liquid ? GC::LRAT
                               : t.mode == Mode::Resv   ? GC::RESV : GC::ORAT;
                ws.production_cmode = Well::ProducerCMode::GRUP;
                ws.group_target.emplace();
                ws.group_target->group_name = system.groups()[t.group].name;
                ws.group_target->target_value = t.value;
                ws.group_target->production_cmode = cmode;
                if (stein_set) {
                    // The tree's own number where it names the same group and mode: those
                    // sum to the group's limit exactly, the inflow lines' only nearly.
                    const auto node = stein_last.find(well.name);
                    if (node != stein_last.end() && node->second.groupTarget.groupName == ws.group_target->group_name
                        && node->second.groupTarget.ctrlMode == cmode && node->second.groupTarget.value > Scalar{0}) {
                        ws.group_target->target_value = node->second.groupTarget.value;
                    }
                }
                holding[ws.group_target->group_name] = cmode;
                ws.group_target_fallback = std::nullopt;
                ws.use_group_target_fallback = false;
                break;
            }
            case Ctrl::Thp:
            case Ctrl::OilRate:
            case Ctrl::Tied: {
                const auto om = own_limit_mode.find(well.name);
                ws.production_cmode = t.control == Ctrl::Thp ? Well::ProducerCMode::THP
                    : (om != own_limit_mode.end() ? om->second : Well::ProducerCMode::ORAT);
                break;
            }
            case Ctrl::Bhp:     ws.production_cmode = Well::ProducerCMode::BHP;  break;
            case Ctrl::Shut:    break;    // handed over below through the dead flag
            default:            break;
            }
            decided.insert(well.name);
            if (ws.production_cmode != cmode_before) {
                switched.insert(well.name);
            }
            assigned_cmode[well.name] = ws.production_cmode;
            dead_now[well.name] = t.control == Ctrl::Shut && !well.shut;
            if (!(well.q_start > Scalar{0}) && rr.well_rate[w] > Scalar{0}) {
                ++this->controller_revivals_[well.name];
            }
            // The other phases where the inflow puts them at the bhp this oil rate needs.
            std::vector<Scalar> q(ws.surface_rates.size(), Scalar{0});
            const Scalar q_oil = rr.well_rate[w];
            if (q_oil > Scalar{0} && well.ipr_b[1] < Scalar{0}) {
                const Scalar bhp = (q_oil - well.ipr_a[1]) / well.ipr_b[1];
                for (int ph = 0; ph < Sys::NP; ++ph) {
                    q[pos[ph]] = -std::max(Sys::ipr(well, ph, bhp), Scalar{0});
                }
            }
            assigned_rates[well.name] = q;
        }
        if (stein_set && !stein_last.empty()) {
            // What was written against the tree it came from: control, holding group, mode,
            // target and rate of every well the tree knows.
            int n = 0, off_control = 0, off_group = 0;
            Scalar off_target = 0, off_rate = 0;
            std::string worst;
            for (int w = 0; w < system.numWells(); ++w) {
                const auto& well = system.wells()[w];
                const auto it = stein_last.find(well.name);
                if (it == stein_last.end() || well.pinned || well.group < 0 || !this->wellState().has(well.name)) {
                    continue;
                }
                ++n;
                const auto& node = it->second;
                const auto& ws = this->wellState().well(well.name);
                const bool theirs = node.modeCategory == ProdNodeModeCategory::Group;
                const bool ours = ws.production_cmode == Well::ProducerCMode::GRUP;
                if (theirs != ours) {
                    ++off_control;
                    worst = well.name;
                }
                const Scalar r = std::abs(rr.well_rate[w] + node.rates[0])
                    / std::max({rr.well_rate[w], -node.rates[0], Scalar{1e-9}});
                if (r > off_rate) { off_rate = r; if (off_control == 0) { worst = well.name; } }
                if (theirs && ours && ws.group_target.has_value()) {
                    const auto& gt = node.groupTarget;
                    off_group += gt.groupName != ws.group_target->group_name
                        || gt.ctrlMode != ws.group_target->production_cmode;
                    off_target = std::max(off_target, std::abs(gt.value - ws.group_target->target_value)
                        / std::max({gt.value, ws.group_target->target_value, Scalar{1e-9}}));
                }
            }
            auto& st = this->controller_stats_;
            ++st.stein_decisions;
            const bool consistent = off_control == 0 && off_group == 0 && off_target <= Scalar{0.01} && off_rate <= Scalar{0.01};
            st.stein_inconsistent += !consistent;
            if (!consistent) {
                deferred_logger.debug(fmt::format(
                    "Controller: written against the balancer's tree under {}: {} wells, {} on another control, "
                    "{} under another group or mode, targets off by {:.1f} %, rates by {:.1f} % ({})", root.name(), n,
                    off_control, off_group, 100.0 * off_target, 100.0 * off_rate, worst));
            }
        }
    }
    // The shut decision handed to the wells: dead while the system holds it shut,
    // and cleared for every other well so a stopped one may be revived.
    for (const auto& name : repaired_dead) { dead_now[name] = true; }
    for (const auto& wp : well_container_) {
        const auto it = dead_now.find(wp->name());
        wp->setNetworkDead(it != dead_now.end() && it->second);
        wp->setNetworkHeld(kept_dead.count(wp->name()) > 0);
    }
    // The group state's controls say what was decided: the holding groups their mode,
    // the groups under one FLD, the rest NONE. Output and legacy's bookkeeping read them.
    for (const auto& g : route_groups) {
        if (!schedule.hasGroup(g, reportStepIdx) || !schedule.getGroup(g, reportStepIdx).isProductionGroup()) {
            continue;
        }
        auto cmode = Group::ProductionCMode::NONE;
        if (const auto it = holding.find(g); it != holding.end()) {
            cmode = it->second;
        } else {
            for (std::string up = g; up != "FIELD" && cmode == Group::ProductionCMode::NONE; ) {
                up = schedule.getGroup(up, reportStepIdx).parent();
                if (holding.count(up)) {
                    cmode = Group::ProductionCMode::FLD;
                }
            }
        }
        this->groupState().production_control(g, cmode);
    }
    // Emit: node pressures and the thp limits they impose, then the wells.
    this->network_.setOwnedNodePressures(new_pressures);
    this->controller_decided_wells_ = decided;
    this->controller_route_wells_ = route_wells;
    controllerMarkDecided_();
    this->controller_assigned_cmode_ = assigned_cmode;
    this->controller_assigned_rates_ = assigned_rates;
    OPM_BEGIN_PARALLEL_TRY_CATCH()
    // As legacy: a well is primed for a control it was switched to, and otherwise left
    // as its solve left it. Priming every time keeps the step from ever converging.
    for (const auto& well : well_container_) {
        if (switched.count(well->name())) {
            well->updateWellStateWithTarget(simulator_, this->groupStateHelper(), this->wellState());
            well->updatePrimaryVariables(this->groupStateHelper());
        }
    }
    OPM_END_PARALLEL_TRY_CATCH("BlackoilWellModel: priming the controller's wells failed: ",
                               simulator_.gridView().comm());
    this->updateAndCommunicateGroupData(reportStepIdx, /*update_wellgrouptarget*/ true);
    this->controller_judge_rejections_ += n_judged_wrong;
    return true;
}

template<typename TypeTag>
std::string
BlackoilWellModel<TypeTag>::
controllerInjectionSignature_() const
{
    std::string sig;
    for (const auto& wp : well_container_) {
        if (!wp->isInjector()) { continue; }
        const auto& ws = this->wellState().well(wp->indexOfWell());
        sig += fmt::format("{}:{}:{};", wp->name(), static_cast<int>(ws.injection_cmode),
                           ws.group_target ? ws.group_target->group_name : std::string{});
    }
    return sig;
}

template<typename TypeTag>
bool
BlackoilWellModel<TypeTag>::
controllerInjectionDecide_(DeferredLogger& deferred_logger, const bool targets_only)
{
    // One Stein tree per injected phase: a group holds the strictest of its injection limits (legacy's
    // target formulas), its injectors share it by guide rate, each capped at its own limits.
    const int step = simulator_.episodeIndex();
    const auto& schedule = this->schedule();
    const auto& summary_state = this->summaryState();
    const auto& pu = this->phaseUsage();
    const Group& field = schedule.getGroup("FIELD", step);
    for (const auto& name : schedule.groupNames(step)) {
        if (schedule.getGroup(name, step).hasSatelliteInjection()) {
            deferred_logger.debug(fmt::format("Controller: injection left to legacy at report step {} ({} has satellite "
                                              "injection)", step, name));
            return false;
        }
    }
    auto& helper = this->groupStateHelper();
    if (param_.group_controller_injection_current_production_) {
        // The current production, not the NUPCOL state the helper holds: after NUPCOL too.
        auto guard = helper.pushWellState(this->wellState());
        helper.updateREINForGroups(field, /*sum_rank=*/this->comm().rank() == 0);
        helper.updateVREPForGroups(field);
        helper.updateReservoirRatesInjectionGroups(field);
        helper.updateSurfaceRatesInjectionGroups(field);
    }
    std::vector<Scalar> group_resv(this->numPhases(), Scalar{0});
    calcInjResvCoeff(/*fipnum=*/0, /*pvtreg=*/0, group_resv);
    using Node = ProdGroupTreeNode<Scalar>;
    using ProdGroupTreeBalancer::Tree;
    struct Phaseinfo { Phase phase; int canonical; InjectorType type; Well::ProducerCMode mode; Group::ProductionCMode gmode; int slot; };
    // Water last: its VREP target subtracts the other phases' injection, decided first.
    const std::array<Phaseinfo, 3> phases{{
        {Phase::GAS, IndexTraits::gasPhaseIdx, InjectorType::GAS, Well::ProducerCMode::GRAT, Group::ProductionCMode::GRAT, 2},
        {Phase::OIL, IndexTraits::oilPhaseIdx, InjectorType::OIL, Well::ProducerCMode::ORAT, Group::ProductionCMode::ORAT, 0},
        {Phase::WATER, IndexTraits::waterPhaseIdx, InjectorType::WATER, Well::ProducerCMode::WRAT, Group::ProductionCMode::WRAT, 1}}};
    for (const auto& wname : schedule.wellNames(step)) {
        const auto& well = schedule.getWell(wname, step);
        if (well.isInjector() && well.injectionControls(summary_state).injector_type == InjectorType::MULTI) {
            return false;
        }
    }
    // What the decision reads about a well, gathered: every rank then decides from the same numbers,
    // and each writes the answer for the wells it owns.
    struct InjWell {
        bool usable{false};
        Well::InjectorCMode cmode{Well::InjectorCMode::CMODE_UNDEFINED};
        std::string holder{};
        Scalar target{0}, bhp{0}, eff{1};
        int pvtreg{0};
        std::vector<Scalar> q, resv, pot;
    };
    const auto& well_names = schedule.wellNames(step);
    const auto& group_names = schedule.groupNames(step);
    const int np = this->numPhases();
    const int record = 7 + 3 * np;
    std::vector<Scalar> buf(well_names.size() * record, Scalar{0});
    std::map<std::string, WellInterface<TypeTag>*> local_wells;
    for (const auto& wp : well_container_) {
        local_wells[wp->name()] = wp.get();
    }
    for (std::size_t i = 0; i < well_names.size(); ++i) {
        const auto& wname = well_names[i];
        const auto index = this->wellState().index(wname);
        if (!index.has_value() || !this->wellState().wellIsOwned(*index, wname)) {
            continue;
        }
        const auto& ws = this->wellState().well(*index);
        const auto it = local_wells.find(wname);
        Scalar* r = buf.data() + i * record;
        r[0] = (ws.status == WellStatus::OPEN && it != local_wells.end() && !it->second->wellIsStopped()) ? 1 : 0;
        r[1] = static_cast<Scalar>(static_cast<int>(ws.injection_cmode));
        if (ws.group_target) {
            const auto g = std::find(group_names.begin(), group_names.end(), ws.group_target->group_name);
            r[2] = g == group_names.end() ? Scalar{0} : static_cast<Scalar>(std::distance(group_names.begin(), g) + 1);
            r[3] = ws.group_target->target_value;
        }
        r[4] = ws.bhp;
        r[5] = ws.efficiency_scaling_factor;
        r[6] = it != local_wells.end() ? static_cast<Scalar>(it->second->pvtRegionIdx()) : Scalar{0};
        for (int p = 0; p < np; ++p) {
            r[7 + p] = ws.surface_rates[p];
            r[7 + np + p] = ws.reservoir_rates[p];
            r[7 + 2 * np + p] = p < static_cast<int>(ws.well_potentials.size()) ? ws.well_potentials[p] : Scalar{0};
        }
    }
    if (this->comm().size() > 1) {
        this->comm().sum(buf.data(), buf.size());
    }
    std::map<std::string, InjWell> gathered;
    for (std::size_t i = 0; i < well_names.size(); ++i) {
        const Scalar* r = buf.data() + i * record;
        InjWell d;
        d.usable = r[0] > Scalar{0.5};
        d.cmode = static_cast<Well::InjectorCMode>(static_cast<int>(std::lround(r[1])));
        const int holder = static_cast<int>(std::lround(r[2]));
        if (holder > 0) {
            d.holder = group_names[holder - 1];
            d.target = r[3];
        }
        d.bhp = r[4];
        d.eff = r[5];
        d.pvtreg = static_cast<int>(std::lround(r[6]));
        d.q.assign(r + 7, r + 7 + np);
        d.resv.assign(r + 7 + np, r + 7 + 2 * np);
        d.pot.assign(r + 7 + 2 * np, r + 7 + 3 * np);
        gathered.emplace(well_names[i], std::move(d));
    }
    std::set<std::string> decided_injectors;
    std::map<std::string, std::map<int, Scalar>> decided_resv;   // group -> active phase -> reservoir rate
    this->controller_injection_moved_ = false;
    for (const auto& ph : phases) {
        if (!pu.phaseIsActive(ph.canonical) || schedule[step].injectionNetwork.has(ph.phase)) {
            continue;       // the injection network stays legacy's for now
        }
        const int pos = pu.canonicalToActivePhaseIdx(ph.canonical);
        Tree<Scalar> tree;
        struct Inj { Well::InjectorCMode cap_mode; };
        std::map<std::string, Inj> injectors;
        std::map<std::string, Group::InjectionCMode> binding;   // group -> the limit it holds
        const auto target = helper.getInjectionGuideTargetMode(ph.phase);
        for (const auto& wname : schedule.wellNames(step)) {
            const auto& well = schedule.getWell(wname, step);
            if (!well.isInjector() || !well.predictionMode()) { continue; }
            const auto controls = well.injectionControls(summary_state);
            if (controls.injector_type != ph.type) { continue; }
            const auto& d = gathered.at(wname);
            if (!d.usable) { continue; }
            Scalar cap = std::numeric_limits<Scalar>::max();
            auto cap_mode = Well::InjectorCMode::BHP;
            if (controls.hasControl(Well::InjectorCMode::RATE) && controls.surface_rate >= 0.0) {
                cap = controls.surface_rate;
                cap_mode = Well::InjectorCMode::RATE;
            }
            if (controls.hasControl(Well::InjectorCMode::RESV) && controls.reservoir_rate >= 0.0) {
                std::vector<Scalar> c(this->numPhases(), Scalar{0});
                calcInjResvCoeff(0, d.pvtreg, c);
                if (c[pos] > Scalar{0} && controls.reservoir_rate / c[pos] < cap) {
                    cap = controls.reservoir_rate / c[pos];
                    cap_mode = Well::InjectorCMode::RESV;
                }
            }
            const Scalar pot = std::abs(d.pot[pos]);
            if (pot > Scalar{0} && pot < cap) {
                cap = pot;
                cap_mode = Well::InjectorCMode::BHP;
            }
            // On its own control a well is what it injects (the tree may still pull it onto GRUP);
            // on GRUP a stale potential below its current rate is no cap.
            const Scalar now = std::abs(d.q[pos]);
            const bool on_grup = d.cmode == Well::InjectorCMode::GRUP && well.isAvailableForGroupControl();
            if (!on_grup || now > cap) {
                cap = now;
                cap_mode = on_grup ? cap_mode : d.cmode;
            }
            // Its own feasibility report: pinned at its bhp limit and short of the share it was given,
            // it cannot do more than it does, whatever its potential says.
            if (param_.group_controller_injection_feedback_
                && on_grup && !d.holder.empty() && d.bhp >= controls.bhp_limit * (1 - Scalar{1e-3})
                && now < d.target * (1 - param_.group_controller_rate_tolerance_)
                && now < cap) {
                cap = now;
                cap_mode = Well::InjectorCMode::BHP;
            }
            Node n;
            n.name = wname;
            n.type = ProdNodeType::Well;
            n.parent = well.groupName();
            n.availableForGroupControl = well.isAvailableForGroupControl();
            n.hasGuideRate = true;
            const auto& g = this->guideRate();
            // Outside group control a well takes no share of its group's target.
            n.fixedGuideRate = !n.availableForGroupControl ? Scalar{0}
                : (g.has(n.name) || g.hasPotentials(n.name))
                ? static_cast<Scalar>(g.get(n.name, target, helper.getWellRateVector(n.name))) : pot;
            if (n.availableForGroupControl && *n.fixedGuideRate <= Scalar{0}) {
                // Just opened, no potential yet: a zero guide would leave its branch without a share.
                n.fixedGuideRate = cap < std::numeric_limits<Scalar>::max() ? cap : Scalar{1};
            }
            n.efficiencyFactor = well.getEfficiencyFactor() * d.eff;
            // At capacity: a group may hold its injectors at zero now.
            const bool grup = on_grup;
            n.rates = {};
            n.rates[ph.slot] = cap < std::numeric_limits<Scalar>::max() ? -cap : -now;
            n.initialRates = n.rates;
            if (cap < std::numeric_limits<Scalar>::max()) { n.Limits[ph.mode] = cap; }
            if (grup) {
                n.mode = Well::ProducerCMode::GRUP;
                n.modeCategory = ProdNodeModeCategory::Group;
            } else {
                n.mode = ph.mode;
                n.modeCategory = ProdNodeModeCategory::Individual;
            }
            tree.emplace(n.name, n);
            injectors.emplace(n.name, Inj{cap_mode});
        }
        if (injectors.empty()) {
            continue;
        }
        // The groups above them, with the strictest of their injection limits for this phase.
        std::function<void(const std::string&)> addGroup = [&](const std::string& gname) {
            if (tree.count(gname)) { return; }
            const auto& group = schedule.getGroup(gname, step);
            Node n;
            n.name = gname;
            n.type = ProdNodeType::Group;
            n.parent = gname == "FIELD" ? std::string{} : group.parent();
            n.availableForGroupControl = gname != "FIELD" && group.injectionGroupControlAvailable(ph.phase);
            n.efficiencyFactor = group.getGroupEfficiencyFactor();
            n.preferredMode = ph.gmode;
            n.hasGuideRate = true;
            if (this->guideRate().has(gname, ph.phase)) {
                n.fixedGuideRate = static_cast<Scalar>(this->guideRate().get(gname, ph.phase));
            }
            Scalar limit = std::numeric_limits<Scalar>::max();
            for (const auto m : {Group::InjectionCMode::RATE, Group::InjectionCMode::RESV, Group::InjectionCMode::REIN,
                                 Group::InjectionCMode::VREP, Group::InjectionCMode::SALE}) {
                // GCONSALE makes SALE the gas mode, as legacy's setCmodeGroup does.
                const bool sale = m == Group::InjectionCMode::SALE && ph.phase == Phase::GAS
                    && schedule[step].gconsale().has(gname);
                const bool own = group.isInjectionGroup() && group.hasInjectionControl(ph.phase)
                    && group.has_control(ph.phase, m) && m != Group::InjectionCMode::SALE;
                if (!(own || sale || group.has_gpmaint_control(ph.phase, m))) {
                    continue;
                }
                Scalar t = helper.injectionGroupTargetForMode(group, ph.phase, group_resv, m);
                if (m == Group::InjectionCMode::RESV || m == Group::InjectionCMode::VREP) {
                    // Its share of the reservoir volume is what the phases already decided leave it.
                    const auto& state_resv = this->groupState().injection_reservoir_rates(gname);
                    Scalar moved = 0;
                    for (const auto& [p_slot, value] : decided_resv[gname]) {
                        if (p_slot != pos) { moved += value - state_resv[p_slot]; }
                    }
                    if (group_resv[pos] > Scalar{0}) { t -= moved / group_resv[pos]; }
                }
                if (t < limit) {
                    limit = std::max(t, Scalar{0});
                    binding[gname] = m;
                }
            }
            if (limit < std::numeric_limits<Scalar>::max()) {
                // A target that follows the current production is damped: after a restart the voidage is
                // transiently far off and the target follows it. GPMAINT and deck rates are not damped.
                const auto m = binding.count(gname) ? binding.at(gname) : Group::InjectionCMode::NONE;
                const bool follows = m == Group::InjectionCMode::REIN || m == Group::InjectionCMode::VREP
                    || m == Group::InjectionCMode::SALE;
                const auto key = fmt::format("{}|{}", gname, static_cast<int>(ph.phase));
                const Scalar d = param_.group_controller_injection_damping_;
                if (const auto it = this->controller_injection_limit_.find(key);
                    follows && d > Scalar{0} && it != this->controller_injection_limit_.end() && it->second > Scalar{0}) {
                    limit = std::clamp(limit, it->second * (1 - d), it->second * (1 + d));
                }
                this->controller_injection_limit_[key] = limit;
                n.Limits[ph.mode] = limit;
            }
            tree.emplace(gname, n);
            if (gname != "FIELD") {
                addGroup(group.parent());
                tree.at(group.parent()).children.push_back(gname);
            }
        };
        for (const auto& [name, inj] : injectors) {
            const auto& gname = tree.at(name).parent;
            addGroup(gname);
            tree.at(gname).children.push_back(name);
        }
        // Injectors outside group control are reductions, as in legacy: their rate comes off every
        // ancestor's limit and they leave the tree (a fixed child in a subgroup would take its share).
        std::map<std::string, Scalar> fixed_rate;
        for (const auto& [name, inj] : injectors) {
            const auto& wn = tree.at(name);
            if (wn.availableForGroupControl) { continue; }
            fixed_rate[name] = -wn.rates[ph.slot];
            Scalar q = -wn.rates[ph.slot] * wn.efficiencyFactor;
            for (std::string g = wn.parent; !g.empty(); g = tree.at(g).parent) {
                auto& gn = tree.at(g);
                if (auto it = gn.Limits.find(ph.mode); it != gn.Limits.end()) {
                    it->second = std::max(it->second - q, Scalar{0});
                }
                q *= gn.efficiencyFactor;
            }
            auto& siblings = tree.at(wn.parent).children;
            siblings.erase(std::remove(siblings.begin(), siblings.end(), name), siblings.end());
            tree.erase(name);
        }
        // Group rates and default guide rates, bottom-up.
        std::function<void(const std::string&)> sum = [&](const std::string& gname) {
            auto& n = tree.at(gname);
            if (n.type != ProdNodeType::Group) { return; }
            std::array<Scalar, 3> r{};
            Scalar guides = 0;
            for (const auto& c : n.children) {
                sum(c);
                const auto& cn = tree.at(c);
                for (int k = 0; k < 3; ++k) { r[k] += cn.efficiencyFactor * cn.rates[k]; }
                guides += cn.efficiencyFactor * cn.fixedGuideRate.value_or(Scalar{0});
            }
            n.rates = r;
            n.initialRates = r;
            if (!n.fixedGuideRate) { n.fixedGuideRate = guides; }
        };
        sum("FIELD");
        DeferredLogger quiet;
        const bool ok = ProdGroupTreeBalancer::balanceTreeForTesting(tree, this->guideRate(),
                                                                     static_cast<Scalar>(param_.group_tree_balancer_tolerance_),
                                                                     std::getenv("OPM_CONTROLLER_TRACE") ? deferred_logger : quiet,
                                                                     /*assignTargets=*/false,
                                                                     /*requireValid=*/false);
        // A held group left over its limit by injectors outside group control is accepted: nothing can be cut.
        if (!ok) {
            deferred_logger.debug(fmt::format("Controller: the injection tree for phase {} is not valid at report step {}; "
                                              "injection left to legacy", static_cast<int>(ph.phase), step));
            return false;
        }
        // Write: held injectors on GRUP at their share under the nearest holding group.
        std::set<std::string> carries;
        std::map<std::string, std::string> holder;
        for (const auto& [name, inj] : injectors) {
            if (!tree.count(name) || tree.at(name).modeCategory != ProdNodeModeCategory::Group) { continue; }
            const auto& n = tree.at(name);
            std::string g = n.parent;
            while (!g.empty()) {
                const auto& gn = tree.at(g);
                if (gn.modeCategory == ProdNodeModeCategory::Individual && gn.Limits.count(ph.mode)) { break; }
                g = gn.parent;
            }
            if (!g.empty()) { holder[name] = g; }
        }
        std::set<std::string> stands;
        if (targets_only) {
            std::map<std::string, std::pair<Scalar, Scalar>> sums;     // group -> (previous, new) held total
            std::set<std::string> fresh;
            for (const auto& [name, g] : holder) {
                const auto& d = gathered.at(name);
                if (d.holder != g) { fresh.insert(g); continue; }
                sums[g].first += d.target;
                sums[g].second += std::max(-tree.at(name).rates[ph.slot], Scalar{0});
            }
            for (const auto& [g, s] : sums) {
                if (!fresh.count(g) && std::abs(s.first - s.second)
                    <= param_.group_controller_rate_tolerance_ * tree.at(g).Limits[ph.mode] + Scalar{1e-7}) {
                    stands.insert(g);
                }
            }
        }
        if (std::getenv("OPM_CONTROLLER_TRACE")) {
            std::string d;
            for (const auto& [name, inj] : injectors) {
                const auto& g = gathered.at(name);
                d += fmt::format(" [{} {} q {:.1f} tgt {:.1f} tree {:.1f} cap {:.1f} holder {} {}]", name,
                                 WellInjectorCMode2String(g.cmode), std::abs(g.q[pos]) * 86400.0,
                                 g.holder.empty() ? -1.0 : g.target * 86400.0,
                                 tree.count(name) ? -tree.at(name).rates[ph.slot] * 86400.0 : -1.0,
                                 tree.count(name) && tree.at(name).Limits.count(ph.mode) ? tree.at(name).Limits.at(ph.mode) * 86400.0 : -1.0,
                                 holder.count(name) ? holder.at(name) : std::string("-"),
                                 holder.count(name) && stands.count(holder.at(name)) ? "stands" : "");
            }
            deferred_logger.debug(fmt::format("CTRLDBG step={} it={} ph{} targets_only={}{}", step,
                                              simulator_.problem().iterationContext().iteration(),
                                              static_cast<int>(ph.phase), targets_only, d));
        }
        std::string trace;
        for (const auto& [name, inj] : injectors) {
            const auto& d = gathered.at(name);
            const auto before = d.cmode;
            // The decision is the same on every rank; only its owner writes it.
            auto* ws = this->wellState().has(name) ? &this->wellState().well(name) : nullptr;
            if (!tree.count(name)) {
                continue;       // outside group control, on its own control
            }
            auto h = holder.find(name);
            const Scalar rate = std::max(-tree.at(name).rates[ph.slot], Scalar{0});
            const Scalar cap = tree.at(name).Limits.count(ph.mode) ? tree.at(name).Limits.at(ph.mode) : rate;
            const bool released = param_.group_controller_injection_feedback_
                && before == Well::InjectorCMode::GRUP && h == holder.end()
                && !d.holder.empty() && tree.count(d.holder);
            if (targets_only && !released
                && (before != Well::InjectorCMode::GRUP || h == holder.end()
                    || d.holder != h->second)) {
                // After NUPCOL the modes stand; only the targets follow production.
                if (before == Well::InjectorCMode::GRUP) { decided_injectors.insert(name); }
                continue;
            }
            const bool at_cap = rate >= cap * (1 - param_.group_controller_rate_tolerance_);
            if (h != holder.end() && before != Well::InjectorCMode::GRUP && at_cap) {
                // Held at its own cap within the tolerance: the switch is not worth an iteration.
                trace += fmt::format(" {} kept {}", name, WellInjectorCMode2String(before));
                continue;
            }
            if (released) {
                // Released at its cap: it stays on GRUP with the cap as its share.
                if (ws != nullptr && ws->group_target) { ws->group_target->target_value = cap; }
                trace += fmt::format(" {} GRUP at cap {:.1f}", name, cap * 86400.0);
                decided_injectors.insert(name);
                continue;
            }
            if (h != holder.end()) {
                const Scalar old = d.holder.empty() ? Scalar{-1} : d.target;
                const bool moved = std::abs(rate - old) > param_.group_controller_rate_tolerance_ * std::abs(rate) + Scalar{1e-7};
                // After NUPCOL a group whose previous shares still meet its target within the tolerance
                // keeps them, or the targets chase Newton's iterate.
                const Scalar share = (targets_only && stands.count(h->second)) ? old : rate;
                this->controller_injection_moved_ = this->controller_injection_moved_ || moved;
                if (ws != nullptr) {
                    ws->injection_cmode = Well::InjectorCMode::GRUP;
                    ws->group_target.emplace();
                    ws->group_target->group_name = h->second;
                    ws->group_target->target_value = share;
                    ws->group_target->injection_cmode = Group::InjectionCMode::NONE;
                    if (const auto b = binding.find(h->second); b != binding.end()) {
                        ws->group_target->injection_cmode = b->second;
                    }
                }
                for (std::string g = tree.at(name).parent; g != h->second; g = tree.at(g).parent) { carries.insert(g); }
                trace += fmt::format(" {} GRUP {:.1f} under {}", name, share * 86400.0, h->second);
            } else if (before == Well::InjectorCMode::GRUP) {
                if (ws != nullptr) {
                    ws->injection_cmode = inj.cap_mode;
                    ws->group_target.reset();
                }
                trace += fmt::format(" {} {} (cap {:.1f}, q {:.1f})", name, WellInjectorCMode2String(inj.cap_mode),
                                     tree.at(name).Limits.count(ph.mode) ? tree.at(name).Limits.at(ph.mode) * 86400.0 : -1.0,
                                     -tree.at(name).rates[ph.slot] * 86400.0);
            }
            decided_injectors.insert(name);
            if (ws != nullptr && ws->injection_cmode != before) {
                auto* w = local_wells.at(name);
                w->updateWellStateWithTarget(simulator_, this->groupStateHelper(), this->wellState());
                w->updatePrimaryVariables(this->groupStateHelper());
            }
        }
        for (const auto& [name, inj] : injectors) {
            if (!this->wellState().has(name)) { continue; }
            const auto& w = this->wellState().well(name);
            trace += fmt::format(" <{}{} {} {:.1f}>", name, tree.count(name) ? "" : " fixed",
                                 WellInjectorCMode2String(w.injection_cmode),
                                 w.group_target ? w.group_target->target_value * 86400.0 : -1.0);
        }
        for (const auto& [gname, gn] : tree) {
            if (gn.type == ProdNodeType::Group && gn.Limits.count(ph.mode)) {
                trace += fmt::format(" [{} limit {:.1f} cat {} q {:.1f}]", gname, gn.Limits.at(ph.mode) * 86400.0,
                                     static_cast<int>(gn.modeCategory), -gn.rates[ph.slot] * 86400.0);
            }
        }
        // The next phase's RESV and VREP see this phase's injection as decided, at each well's own
        // reservoir/surface ratio (a field-average gas factor was 22 % off). Held here, not in the group
        // state: those containers hold each rank's own share until communicate_rates sums them.
        std::map<std::string, Scalar> resv_sum;
        for (const auto& [name, inj] : injectors) {
            const auto& d = gathered.at(name);
            const auto& well = schedule.getWell(name, step);
            const Scalar q = tree.count(name) ? -tree.at(name).rates[ph.slot] : fixed_rate[name];
            const Scalar s = std::abs(d.q[pos]);
            const Scalar ratio = s > Scalar{1e-12} ? std::abs(d.resv[pos]) / s : group_resv[pos];
            Scalar r = q * ratio * well.getEfficiencyFactor() * d.eff;
            for (std::string g = well.groupName(); tree.count(g); g = tree.at(g).parent) {
                resv_sum[g] += r;
                r *= tree.at(g).efficiencyFactor;
            }
        }
        for (const auto& [gname, r] : resv_sum) {
            decided_resv[gname][pos] = r;
        }
        for (const auto& [well, g] : holder) {
            if (const auto b = binding.find(g); b != binding.end()) {
                this->groupState().injection_control(g, ph.phase, b->second);
            }
        }
        for (const auto& g : carries) {
            this->groupState().injection_control(g, ph.phase, Group::InjectionCMode::FLD);
        }
        if (std::getenv("OPM_CONTROLLER_TRACE") != nullptr && !trace.empty()) {
            deferred_logger.debug(fmt::format("CTRLTRACE step={} it={} injection phase {}:{}", step,
                                              simulator_.problem().iterationContext().iteration(),
                                              static_cast<int>(ph.phase), trace));
        }
    }
    this->controller_injection_decided_ = decided_injectors;
    controllerMarkDecided_();
    return true;
}

} // namespace Opm

#endif // OPM_BLACKOILWELLMODEL_CONTROLLER_IMPL_HEADER_INCLUDED
