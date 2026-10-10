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

#ifndef OPM_BLACKOILWELLMODEL_FACILITY_REFERENCE_IMPL_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_FACILITY_REFERENCE_IMPL_HEADER_INCLUDED

// Which open/shut states of the network's transition wells are facility solutions at this
// reservoir state, by the wells' own equations and the network alone. Every combination of the
// candidates is driven to a fixed point of (rates at the node pressures, node pressures at the
// rates); the state handed over is scored against the set found. Diagnostic only, serial only;
// the run is left exactly as it was.

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace Opm {

template<typename TypeTag>
void
BlackoilWellModel<TypeTag>::
facilityReference_(DeferredLogger& deferred_logger)
{
    const int step = simulator_.episodeIndex();
    const auto& schedule = this->schedule();
    const auto& network = schedule[step].network();
    if (this->comm().size() > 1) {
        return;
    }
    // Without a network every thp well's pressure is its own limit; the states are still coupled
    // through the group limits and lift gas as handed over.
    const bool with_network = network.active() && this->network_.active();
    const auto& pu = this->phaseUsage();
    if (!pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
        return;
    }
    static const int max_cand = [] {
        const char* e = std::getenv("OPM_FACILITY_REFERENCE_MAX");
        return e ? std::max(1, std::atoi(e)) : 5;
    }();
    static const auto step_range = []() -> std::pair<int, int> {
        const char* e = std::getenv("OPM_FACILITY_REFERENCE_STEPS");
        if (!e) { return {0, 1 << 30}; }
        int a = 0, b = 1 << 30;
        std::sscanf(e, "%d-%d", &a, &b);
        return {a, b};
    }();
    static const bool detail = std::getenv("OPM_FACILITY_REFERENCE_DETAIL") != nullptr;
    // OPM_FACILITY_REFERENCE_PROBE=1: stopped wells asked on a probe (the facility check's way) instead of their
    // own object; the two can disagree, and the controller's runs are better with the object.
    static const bool use_probe = std::getenv("OPM_FACILITY_REFERENCE_PROBE") != nullptr;
    // OPM_FACILITY_REFERENCE_FRACTIONS=inherited|explicit|implicit: the lift test's tubing fractions.
    static const LiftFractions fractions = liftFractions_(std::getenv("OPM_FACILITY_REFERENCE_FRACTIONS")
                                                          ? std::getenv("OPM_FACILITY_REFERENCE_FRACTIONS") : "inherited");
    if (step < step_range.first || step > step_range.second) {
        return;
    }
    constexpr int max_fixed_point = 14;
    const int oil = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
    const Scalar tol_pressure = param_.group_controller_network_tolerance_;
    const auto& iterCtx = simulator_.problem().iterationContext();
    const auto& summary_state = this->summaryState();
    using ThpPoint = typename WellInterface<TypeTag>::ThpPoint;

    struct Cand
    {
        std::string name, node;
        Scalar own_thp{0};
        bool flowing{false}, flowing_prev{false}, live_stopped{false}, closed{false};
        const WellInterface<TypeTag>* live{nullptr};
        int rank() const
        {
            return flowing != flowing_prev ? 0 : !live_stopped && !closed ? 1 : live_stopped ? 2 : 3;
        }
    };
    std::vector<Cand> cands;
    std::vector<Cand> followers;   // flowing on thp, not candidates: re-solved where their node moves
    int deck_kept = 0, closed_by_run = 0;
    for (const auto& name : schedule.wellNames(step)) {
        const auto& w = schedule.getWell(name, step);
        if (!w.isProducer() || w.getStatus() != Well::Status::OPEN || !w.predictionMode()
            || w.vfp_table_number() <= 0 || !this->wellState().has(name)) {
            continue;
        }
        const auto controls = w.productionControls(summary_state);
        if (with_network ? !network.has_node(w.groupName()) : !controls.hasControl(Well::ProducerCMode::THP)) {
            continue;
        }
        const auto& ws = this->wellState().well(name);
        const auto live = std::find_if(well_container_.begin(), well_container_.end(),
                                       [&name](const auto& wp) { return wp->name() == name; });
        Cand c;
        c.name = name;
        c.node = w.groupName();
        c.own_thp = with_network ? Scalar(0) : controls.thp_limit;
        c.live = live != well_container_.end() ? live->get() : nullptr;
        c.live_stopped = c.live != nullptr && c.live->wellIsStopped();
        c.closed = c.live == nullptr && ws.status == Well::Status::SHUT;
        if (c.live == nullptr) {
            // Closed by the run (WTEST, economic limit, "cannot operate" at a step's end): no decider but the
            // deck can bring it back, so it is not this step's decision.
            closed_by_run += c.closed;
            continue;
        }
        c.flowing = c.live != nullptr && !c.live_stopped && -ws.surface_rates[oil] > Scalar(0);
        c.flowing_prev = c.flowing;
        if (this->prevWellState().has(name)) {
            const auto& pws = this->prevWellState().well(name);
            c.flowing_prev = pws.status == Well::Status::OPEN && -pws.surface_rates[oil] > Scalar(0);
        }
        // Shut before this step: WTEST decides when it is tried again, not the facility.
        const bool prev_open = !this->prevWellState().has(name)
            || this->prevWellState().well(name).status == Well::Status::OPEN;
        if (!c.flowing && !prev_open) {
            ++deck_kept;
            continue;
        }
        if (!c.flowing || c.flowing != c.flowing_prev) {
            cands.push_back(c);
        } else if (with_network && ws.production_cmode == Well::ProducerCMode::THP) {
            followers.push_back(c);
        }
    }
    if (cands.empty()) {
        deferred_logger.debug(fmt::format("Facility reference: step {} iteration {}: no candidates (deck-kept {} closed-by-run {})",
                                          step, iterCtx.iteration(), deck_kept, closed_by_run));
        return;
    }
    std::stable_sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return a.rank() < b.rank(); });
    const int dropped = std::max(0, static_cast<int>(cands.size()) - max_cand);
    if (dropped > 0) {
        cands.resize(max_cand);
    }
    const int n = static_cast<int>(cands.size());

    // From here the state is written to; everything is put back below.
    const auto saved_ws = this->wellState();
    const auto saved_active = this->active_wgstate_;
    const auto saved_nupcol = this->nupcol_wgstate_;
    const auto hand_over_p = with_network ? this->network_.nodePressures() : std::map<std::string, Scalar>{};

    struct Combo
    {
        unsigned mask{0};
        bool converged{false}, dead_open{false}, undetermined{false}, undetermined_shut{false}, cascade{false};
        bool group_over{false}, held_ok{false};
        std::string held_note;
        std::vector<std::string> excused;   // shut wells the group could neither take free nor hold
        int iterations{0};
        Scalar oil{0}, residual{0};
        std::vector<std::string> could_flow;
        std::map<std::string, Scalar> p;
        bool strict() const
        {
            return converged && could_flow.size() == excused.size() && !undetermined_shut && (!group_over || held_ok);
        }
    };
    std::vector<Combo> combos;
    std::string failure;
    int solves_total = 0, scans_used = 0, scan_solves = 0, unstable_fixed = 0;
    const int np = pu.numActivePhases();
    {
        auto log_guard = this->groupStateHelper().pushLogger(/*do_mpi_gather*/ false);
        DeferredLogger scratch;
        try {
            // Stopped candidates are asked on a probe (the lift test does that itself; the held test here too).
            std::map<std::string, WellInterfacePtr> probes;
            for (const auto& c : cands) {
                if (use_probe && c.live_stopped) {
                    probes[c.name] = this->makeProbe_(c.name, scratch);
                }
            }
            auto wellOf = [&](const Cand& c) -> const WellInterface<TypeTag>* {
                return use_probe && c.live_stopped ? probes.at(c.name).get() : c.live;
            };
            auto ratesAt = [&](const WellInterface<TypeTag>* wi, const Scalar bhp) {
                ++scan_solves;
                return this->inflowAtBhp_(*wi, bhp, scratch);
            };
            // The one lift test (liftTest_): the well's equations under thp control on the stable branch, the
            // scan where they do not settle.
            auto solveAt = [&](const Cand& c, const Scalar thp) -> ThpPoint {
                const WellInterface<TypeTag>* wi = wellOf(c);
                const auto a = this->liftTest_(*wi, thp, wi->getALQ(this->wellState()), scratch, use_probe, fractions);
                solves_total += a.solves;
                scans_used += a.by_scan;
                unstable_fixed += a.unstable;
                ThpPoint pt;
                pt.thp = thp;
                pt.determined = a.determined;
                pt.lifts = a.lifts && a.stable;
                pt.bhp = a.bhp;
                pt.flux = a.flux;
                pt.rates = a.rates;
                return pt;
            };
            auto nodeP = [&](const std::map<std::string, Scalar>& p, const Cand& c) {
                if (!with_network) { return c.own_thp; }
                const auto it = p.find(c.node);
                return it != p.end() ? it->second : Scalar(0);
            };

            // The worst production group over a RATE-action limit at the rates in the well state.
            struct GroupOver { bool over{false}; std::string group; int phase{-1}; Scalar target{0}, current{0}; };
            auto worstGroupOver = [&]() {
                GroupOver worst;
                const Scalar tol = param_.group_controller_rate_tolerance_;
                const auto& helper = this->groupStateHelper();
                for (const auto& gname : schedule.groupNames(step)) {
                    const auto& group = schedule.getGroup(gname, step);
                    if (!group.isProductionGroup()) { continue; }
                    const auto controls = group.productionControls(summary_state);
                    const auto& action = controls.group_limit_action;
                    const bool all_rate = action.allRates == Group::ExceedAction::RATE;
                    auto live = [&](const int canonical) {
                        if (!pu.phaseIsActive(canonical)) { return Scalar(0); }
                        return helper.sumWellSurfaceRates(group, pu.canonicalToActivePhaseIdx(canonical), false);
                    };
                    const Scalar o = live(IndexTraits::oilPhaseIdx), w = live(IndexTraits::waterPhaseIdx), g = live(IndexTraits::gasPhaseIdx);
                    auto consider = [&](const Group::ProductionCMode mode, const Group::ExceedAction act, const Scalar target,
                                        const Scalar current, const int phase) {
                        if (!group.has_control(mode) || !(all_rate || act == Group::ExceedAction::RATE) || !(target > Scalar(0))) { return; }
                        if (current > (Scalar(1) + tol) * target && current / target > (worst.over ? worst.current / worst.target : Scalar(0))) {
                            worst = {true, gname, phase, target, current};
                        }
                    };
                    consider(Group::ProductionCMode::ORAT, action.oil, controls.oil_target, o, IndexTraits::oilPhaseIdx);
                    consider(Group::ProductionCMode::WRAT, action.water, controls.water_target, w, IndexTraits::waterPhaseIdx);
                    consider(Group::ProductionCMode::GRAT, action.gas, controls.gas_target, g, IndexTraits::gasPhaseIdx);
                    consider(Group::ProductionCMode::LRAT, action.liquid, controls.liquid_target, o + w, -1);
                }
                return worst;
            };
            auto underGroup = [&](const std::string& well, const std::string& group) {
                std::string g = schedule.getWell(well, step).groupName();
                while (true) {
                    if (g == group) { return true; }
                    if (g == "FIELD" || g.empty()) { return false; }
                    g = schedule.getGroup(g, step).parent();
                }
            };
            auto phaseRate = [&](const std::vector<Scalar>& r, const int phase) {
                if (phase < 0) {
                    return -(r[oil] + (pu.phaseIsActive(IndexTraits::waterPhaseIdx) ? r[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)] : Scalar(0)));
                }
                return pu.phaseIsActive(phase) ? -r[pu.canonicalToActivePhaseIdx(phase)] : Scalar(0);
            };
            // Over a group limit: the excess taken off the open candidates under that group in proportion. Each
            // is feasible held where its tubing lifts the reduced rate at its thp (margin >= 0 at the bhp that
            // gives the rate). Fixed-thp wells only; with a network the node pressures would move too.
            auto heldState = [&](Combo& cb, const unsigned mask) {
                const auto worst = worstGroupOver();
                if (!worst.over) { return; }
                std::vector<int> under;
                Scalar cand_sum = 0;
                for (int i = 0; i < n; ++i) {
                    if ((mask & (1u << i)) != 0 && underGroup(cands[i].name, worst.group)) {
                        under.push_back(i);
                        cand_sum += phaseRate(this->wellState().well(cands[i].name).surface_rates, worst.phase);
                    }
                }
                const Scalar excess = worst.current - worst.target;
                if (under.empty() || !(cand_sum > excess)) {
                    cb.held_note = fmt::format(" held-infeasible({}: excess {:.0f} > candidates {:.0f})", worst.group,
                                               excess * 86400.0, cand_sum * 86400.0);
                    return;
                }
                const Scalar f = Scalar(1) - excess / cand_sum;
                Scalar oil_off = 0;
                for (const int i : under) {
                    const auto& c = cands[i];
                    const WellInterface<TypeTag>* wi = wellOf(c);
                    const Scalar alq = wi->getALQ(this->wellState());
                    const auto& ws = this->wellState().well(c.name);
                    const Scalar q_oil = -ws.surface_rates[oil] * f;
                    // bhp giving the reduced oil rate: rates fall with bhp, so bisect from the operating bhp upward
                    Scalar a = ws.bhp, b = ws.bhp + 10.0e5;
                    std::vector<Scalar> r;
                    for (int k = 0; k < 8; ++k) {
                        r = ratesAt(wi, b);
                        if (-r[oil] <= q_oil) { break; }
                        a = b;
                        b += (b - ws.bhp);
                    }
                    for (int k = 0; k < 12 && b - a > 0.01e5; ++k) {
                        const Scalar mid = 0.5 * (a + b);
                        r = ratesAt(wi, mid);
                        if (-r[oil] > q_oil) { a = mid; } else { b = mid; }
                    }
                    const Scalar m = this->thpMarginAt_(*wi, b, nodeP(cb.p, c), alq, r, scratch);
                    ++scan_solves;
                    if (!(m >= Scalar(0))) {
                        cb.held_note = fmt::format(" held-infeasible({} at {:.0f} sm3/d: margin {:.2f} bar; bhp {:.2f} from {:.2f}, "
                                                   "rates w/o/g {:.0f}/{:.0f}/{:.0f}, well thp {:.2f}, alq {:.0f})",
                                                   c.name, q_oil * 86400.0, m * 1e-5, b * 1e-5, ws.bhp * 1e-5,
                                                   -r[0] * 86400.0, -r[1] * 86400.0, -r[2] * 86400.0, ws.thp * 1e-5, alq * 86400.0);
                        return;
                    }
                    oil_off += -ws.surface_rates[oil] * (Scalar(1) - f);
                }
                cb.held_ok = true;
                cb.oil -= oil_off;
                cb.held_note = fmt::format(" held({} x{:.2f})", worst.group, f);
            };

            for (unsigned mask = 0; mask < (1u << n); ++mask) {
                Combo cb;
                cb.mask = mask;
                auto p = hand_over_p;
                std::map<std::string, Scalar> follower_at;   // node pressure each follower was last solved at
                for (int k = 0; k < max_fixed_point; ++k) {
                    ++cb.iterations;
                    for (int i = 0; i < n && !cb.dead_open && !cb.undetermined; ++i) {
                        const auto& c = cands[i];
                        auto& ws = this->wellState().well(c.name);
                        if ((mask & (1u << i)) == 0) {
                            std::fill(ws.surface_rates.begin(), ws.surface_rates.end(), Scalar(0));
                            continue;
                        }
                        const auto pt = solveAt(c, nodeP(p, c));
                        if (!pt.determined) {
                            cb.undetermined = true;
                        } else if (!pt.lifts) {
                            cb.dead_open = true;
                        } else {
                            ws.open();
                            ws.surface_rates = pt.flux;
                        }
                    }
                    if (cb.dead_open || cb.undetermined) {
                        break;
                    }
                    for (const auto& f : followers) {
                        const Scalar pf = nodeP(p, f);
                        const Scalar last = follower_at.count(f.name) ? follower_at.at(f.name) : nodeP(hand_over_p, f);
                        if (std::abs(pf - last) <= tol_pressure) {
                            continue;
                        }
                        const auto pt = solveAt(f, pf);
                        follower_at[f.name] = pf;
                        if (!pt.determined) {
                            continue;
                        }
                        if (!pt.lifts) {
                            cb.cascade = true;
                            break;
                        }
                        this->wellState().well(f.name).surface_rates = pt.flux;
                    }
                    if (cb.cascade) {
                        break;
                    }
                    if (!with_network) {
                        cb.converged = true;
                        break;
                    }
                    this->updateAndCommunicateGroupData(step, /*update_wellgrouptarget*/ false);
                    const auto p_new = this->network_.pressuresAtGroupRates(step);
                    // Damped after the first steps: rates and node pressures alternate on a steep tubing curve.
                    const Scalar omega = k < 2 ? Scalar(1) : Scalar(0.5);
                    Scalar diff = 0;
                    for (const auto& [node, pn] : p_new) {
                        if (const auto it = p.find(node); it != p.end()) {
                            diff = std::max(diff, std::abs(pn - it->second));
                            it->second += omega * (pn - it->second);
                        } else {
                            p[node] = pn;
                        }
                    }
                    cb.residual = diff;
                    if (diff <= tol_pressure) {
                        cb.converged = true;
                        break;
                    }
                }
                if (cb.converged) {
                    cb.group_over = this->controllerGroupLimitViolated_();
                    for (const auto& name : schedule.wellNames(step)) {
                        if (this->wellState().has(name) && schedule.getWell(name, step).isProducer()) {
                            cb.oil += std::max(-this->wellState().well(name).surface_rates[oil], Scalar(0));
                        }
                    }
                    cb.p = p;
                    if (cb.group_over && !with_network) {
                        heldState(cb, mask);
                    }
                    for (int i = 0; i < n; ++i) {
                        if ((mask & (1u << i)) != 0) {
                            continue;
                        }
                        const auto pt = solveAt(cands[i], nodeP(p, cands[i]));
                        if (pt.determined && pt.lifts) {
                            cb.could_flow.push_back(fmt::format("{}(alq {:.0f})", cands[i].name,
                                                                wellOf(cands[i])->getALQ(this->wellState()) * 86400.0));
                        } else if (!pt.determined) {
                            cb.undetermined_shut = true;
                        }
                    }
                }
                combos.push_back(std::move(cb));
                this->wellState() = saved_ws;
            }
            for (auto& cb : combos) {
                if (!cb.converged) { continue; }
                for (int i = 0; i < n; ++i) {
                    if ((cb.mask & (1u << i)) != 0) { continue; }
                    const bool listed = std::any_of(cb.could_flow.begin(), cb.could_flow.end(),
                        [&](const std::string& f) { return f.rfind(cands[i].name + "(", 0) == 0; });
                    const auto& opened = combos[cb.mask | (1u << i)];
                    if (listed && opened.converged && opened.group_over && !opened.held_ok) {
                        cb.excused.push_back(cands[i].name);
                    }
                }
            }
        } catch (const std::exception& e) {
            failure = e.what();
        }
        log_guard.discard();
    }
    this->wellState() = saved_ws;
    this->active_wgstate_ = saved_active;
    this->nupcol_wgstate_ = saved_nupcol;

    unsigned chosen = 0, prev_mask = 0;
    for (int i = 0; i < n; ++i) {
        chosen |= cands[i].flowing ? (1u << i) : 0u;
        prev_mask |= cands[i].flowing_prev ? (1u << i) : 0u;
    }
    auto bits = [&](const unsigned m) {
        std::string s;
        for (int i = 0; i < n; ++i) {
            s += (i ? " " : "") + cands[i].name + ((m & (1u << i)) ? "+" : "-");
        }
        return s;
    };
    auto popcount = [](unsigned m) { int c = 0; for (; m; m &= m - 1) { ++c; } return c; };
    const Combo* least = nullptr;   // least change from the step's start, then most oil
    const Combo* most = nullptr;    // most oil
    int n_strict = 0, n_adm = 0, n_undet = 0, n_undet_shut = 0, n_nofp = 0, n_cascade = 0, n_dead_open = 0, n_group_over = 0;
    for (const auto& cb : combos) {
        n_adm += cb.converged;
        n_group_over += cb.converged && cb.group_over;
        n_undet += cb.undetermined;
        n_undet_shut += cb.converged && cb.undetermined_shut;
        n_cascade += cb.cascade;
        n_dead_open += cb.dead_open;
        n_nofp += !cb.converged && !cb.dead_open && !cb.undetermined && !cb.cascade;
        if (!cb.strict()) {
            continue;
        }
        ++n_strict;
        if (!most || cb.oil > most->oil) {
            most = &cb;
        }
        if (!least || popcount(cb.mask ^ prev_mask) < popcount(least->mask ^ prev_mask)
            || (popcount(cb.mask ^ prev_mask) == popcount(least->mask ^ prev_mask) && cb.oil > least->oil)) {
            least = &cb;
        }
    }
    const Combo* chosen_cb = chosen < combos.size() ? &combos[chosen] : nullptr;
    const bool chosen_strict = chosen_cb && chosen_cb->strict();
    // A cliff: the chosen state is consistent apart from shut wells that lift at its node pressures,
    // and every one of them dies once opened (the node pressure rises with its own flow).
    bool cliff = n_strict == 0 && chosen_cb && chosen_cb->converged && !chosen_cb->could_flow.empty();
    for (int i = 0; cliff && i < n; ++i) {
        if ((chosen & (1u << i)) == 0 && std::any_of(chosen_cb->could_flow.begin(), chosen_cb->could_flow.end(),
                [&](const std::string& f) { return f.rfind(cands[i].name + "(", 0) == 0; })) {
            const auto& opened = combos[chosen | (1u << i)];
            cliff = cliff && (opened.dead_open || !opened.converged);
        }
    }
    const std::string verdict = !failure.empty() ? "FAILED"
        : cliff ? "CLIFF"
        : n_strict == 0 ? (chosen_cb && chosen_cb->converged ? "NONE-STRICT" : "NONE")
        : n_strict == 1 ? (chosen_strict ? "UNIQUE-OK" : "UNIQUE-MISSED")
        : chosen_strict ? (chosen_cb == least ? "MULTI-LEAST" : chosen_cb == most ? "MULTI-MOST" : "MULTI-OTHER")
        : "MULTI-MISSED";
    auto describe = [&](const Combo* cb) {
        if (!cb) { return std::string("-"); }
        return fmt::format("[{}] oil {:.0f}{}{}", bits(cb->mask), cb->oil * 86400.0,
                           cb->converged ? (std::string(cb->group_over ? " group-over" : "") + cb->held_note
                                            + (cb->undetermined_shut ? " shut-undetermined" : "")
                                            + (cb->excused.empty() ? "" : " excused: " + fmt::format("{}", fmt::join(cb->excused, ","))))
                                         : cb->dead_open ? " dead-open" : cb->undetermined ? " undetermined"
                                         : cb->cascade ? " cascade" : fmt::format(" no-fixed-point ({:.2f} bar)", cb->residual * 1e-5),
                           cb->could_flow.empty() ? "" : " could flow: " + fmt::format("{}", fmt::join(cb->could_flow, ",")));
    };
    deferred_logger.debug(fmt::format(
        "Facility reference: step {} iteration {}: candidates {} [{}]{} followers {} deck-kept {} closed-by-run {} | combos {} admissible {} strict {} "
        "| verdict {} | chosen {} | least-change {} | most-oil {} | dead-open {} undetermined-open {} "
        "undetermined-shut {} no-fixed-point {} cascade {} group-over {} | well solves {} scans {} unstable {} ({} solves){}",
        step, iterCtx.iteration(), n, bits(chosen), dropped ? fmt::format(" (+{} dropped)", dropped) : "",
        followers.size(), deck_kept, closed_by_run, combos.size(), n_adm, n_strict, verdict, describe(chosen_cb), describe(least), describe(most),
        n_dead_open, n_undet, n_undet_shut, n_nofp, n_cascade, n_group_over, solves_total, scans_used, unstable_fixed, scan_solves,
        failure.empty() ? "" : " | failed: " + failure));
    if (detail) {
        for (const auto& cb : combos) {
            deferred_logger.debug(fmt::format("Facility reference detail: step {} iteration {}: {} iterations {}{}",
                                              step, iterCtx.iteration(), describe(&cb), cb.iterations,
                                              cb.strict() ? " STRICT" : ""));
        }
    }
}

} // namespace Opm

#endif
