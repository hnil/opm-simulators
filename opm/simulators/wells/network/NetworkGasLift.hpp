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
#ifndef OPM_NETWORK_GAS_LIFT_HEADER_INCLUDED
#define OPM_NETWORK_GAS_LIFT_HEADER_INCLUDED

#include <opm/simulators/wells/network/NetworkReducedSolve.hpp>
#include <opm/simulators/wells/network/NetworkSolve.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <istream>
#include <limits>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace Opm::NetworkSolve {

/// The system's oil and produced gas, and their derivatives with respect to each lifted well's lift gas,
/// at a converged answer of the reduced route.
template<class Scalar>
struct LiftGradients
{
    bool ok = false;
    Scalar oil = 0, gas = 0;
    std::vector<Scalar> doil, dgas;   // per entry of the wells asked about, per unit lift gas
};

/// dF/da_w = dF/da_w|_p - lambda^T dR/da_w with J^T lambda = dF/dp: the node pressures' response through
/// the node rows R(p, a) = 0, so a well's gradient includes what its lift gas does to its node and to every
/// other well, the group allocation included (it is inside R). Partial derivatives by differences of the
/// reduced residual: nodes + wells + 1 evaluations, no solve.
template<class Sys>
LiftGradients<typename Sys::ScalarType>
liftGradients(Sys& system, const std::vector<typename Sys::ScalarType>& p, const std::vector<int>& wells,
              const typename Sys::ScalarType h_alq, const typename Sys::ScalarType h_p)
{
    using Scalar = typename Sys::ScalarType;
    LiftGradients<Scalar> out;
    const int n = system.numNodes();
    auto totals = [&] {
        const auto& x = system.reducedState();
        Scalar oil = 0, gas = 0;
        for (int w = 0; w < system.numWells(); ++w) {
            oil += system.wells()[w].efficiency * x[system.qwIdx(w, 1)];
            gas += system.wells()[w].efficiency * x[system.qwIdx(w, 2)];
        }
        return std::pair{oil, gas};
    };
    const auto r0 = system.reducedResidual(p);
    const auto [f0, g0] = totals();
    out.oil = f0;
    out.gas = g0;
    DenseMatrix<Scalar> jt(n);
    std::vector<Scalar> dfdp(n), dgdp(n);
    for (int j = 0; j < n; ++j) {
        auto pj = p;
        pj[j + 1] += h_p;
        const auto rj = system.reducedResidual(pj);
        const auto [fj, gj] = totals();
        for (int i = 0; i < n; ++i) { jt(j, i) = (rj[i] - r0[i]) / h_p; }
        dfdp[j] = (fj - f0) / h_p;
        dgdp[j] = (gj - g0) / h_p;
    }
    std::vector<Scalar> lf, lg;
    if (n > 0 && (!jt.solve(dfdp, lf) || !jt.solve(dgdp, lg))) {
        (void)system.reducedResidual(p);
        return out;
    }
    for (const int w : wells) {
        const Scalar a0 = system.wells()[w].alq;
        system.setWellAlq(w, a0 + h_alq);
        const auto rw = system.reducedResidual(p);
        const auto [fw, gw] = totals();
        system.setWellAlq(w, a0);
        Scalar df = (fw - f0) / h_alq, dg = (gw - g0) / h_alq;
        for (int i = 0; i < n; ++i) {
            const Scalar dr = (rw[i] - r0[i]) / h_alq;
            df -= lf[i] * dr;
            dg -= lg[i] * dr;
        }
        out.doil.push_back(df);
        out.dgas.push_back(dg);
    }
    (void)system.reducedResidual(p);
    out.ok = true;
    return out;
}

/// One lifted well: its index in the system, its lift-gas bounds and LIFTOPT weights.
template<class Scalar>
struct LiftWell { int w; Scalar lo, hi, wf, gf; };

/// A group's lift-gas limit (GLIFTOPT), net of satellite lift gas, over the lifted wells below it.
template<class Scalar>
struct LiftCap { Scalar cap; std::vector<int> wells; };

template<class Scalar>
struct LiftProblem
{
    Scalar inc = 0, eco = 0;
    /// A lifted well dead when the allocation starts comes back only through lift gas: if it ends with none, it is
    /// shut again (else the route's own lines may keep it flowing where the well model could not lift it).
    bool revive_only_by_lift = true;
    std::vector<LiftWell<Scalar>> wells;
    std::vector<LiftCap<Scalar>> caps;
};

enum class LiftMethod { Trials, Steps, Projected };

inline LiftMethod liftMethodFromString(const std::string& s)
{
    if (s == "trials") { return LiftMethod::Trials; }
    if (s == "projected") { return LiftMethod::Projected; }
    return LiftMethod::Steps;
}

struct LiftStats { int solves = 0, gradients = 0, moves = 0; };

/// The economic objective the methods are measured by: the system's oil less LIFTOPT's minimum gradient times
/// the lift gas -- an allocation is optimal when no move improves it within the limits.
template<class Sys, class Result>
typename Sys::ScalarType liftObjective(const Sys& system, const Result& r, const LiftProblem<typename Sys::ScalarType>& prob)
{
    using Scalar = typename Sys::ScalarType;
    Scalar oil = 0, lift = 0;
    for (int k = 0; k < system.numWells(); ++k) { oil += system.wells()[k].efficiency * r.well_rate[k]; }
    for (const auto& l : prob.wells) { lift += system.wells()[l.w].alq; }
    return oil - prob.eco * lift;
}

template<class Sys>
bool withinCaps(const Sys& system, const LiftProblem<typename Sys::ScalarType>& prob)
{
    for (const auto& c : prob.caps) {
        typename Sys::ScalarType total = 0;
        for (const int w : c.wells) { total += system.wells()[w].alq; }
        if (total > c.cap * (1 + 1e-9)) { return false; }
    }
    return true;
}

/// Lift gas onto the wells' bounds and the groups' limits (water-filling: one shift per limit, clamped).
template<class Scalar>
void projectLift(std::vector<Scalar>& a, const LiftProblem<Scalar>& prob)
{
    for (std::size_t k = 0; k < a.size(); ++k) { a[k] = std::clamp(a[k], prob.wells[k].lo, prob.wells[k].hi); }
    for (const auto& c : prob.caps) {
        std::vector<std::size_t> in;
        for (std::size_t k = 0; k < prob.wells.size(); ++k) {
            if (std::find(c.wells.begin(), c.wells.end(), prob.wells[k].w) != c.wells.end()) { in.push_back(k); }
        }
        auto sum = [&](const Scalar tau) {
            Scalar t = 0;
            for (const auto k : in) { t += std::clamp(a[k] - tau, prob.wells[k].lo, prob.wells[k].hi); }
            return t;
        };
        if (sum(0) <= c.cap) { continue; }
        Scalar lo = 0, hi = 0;
        for (const auto k : in) { hi = std::max(hi, a[k] - prob.wells[k].lo); }
        for (int it = 0; it < 60; ++it) {
            const Scalar mid = (lo + hi) / 2;
            (sum(mid) > c.cap ? lo : hi) = mid;
        }
        for (const auto k : in) { a[k] = std::clamp(a[k] - hi, prob.wells[k].lo, prob.wells[k].hi); }
    }
}

/// The lift-gas allocation of one decision on the route's converged answer r, on the frozen IPRs. solve(start,
/// keep_dead) re-solves the route; held(name) says a well is not to move this step. Every method first gives a
/// dead well the least lift gas that makes it flow, if its oil per gas pays.
template<class Sys, class Result, class Solve, class Held>
LiftStats allocateLiftGasCore(Sys& system, Result& r, const LiftProblem<typename Sys::ScalarType>& prob,
                              const LiftMethod method, Solve&& solve, Held&& held);

template<class Sys, class Result, class Solve, class Held>
LiftStats allocateLiftGas(Sys& system, Result& r, const LiftProblem<typename Sys::ScalarType>& prob,
                          const LiftMethod method, Solve&& solve, Held&& held)
{
    using Scalar = typename Sys::ScalarType;
    std::vector<int> dead_at_start;
    for (const auto& l : prob.wells) { if (!(r.well_rate[l.w] > Scalar{0})) { dead_at_start.push_back(l.w); } }
    auto st = allocateLiftGasCore(system, r, prob, method, solve, held);
    if (prob.revive_only_by_lift && r.converged) {
        bool shut = false;
        for (const int w : dead_at_start) {
            if (r.well_rate[w] > Scalar{0} && !(system.wells()[w].alq > Scalar{0})) { system.killWell(w); shut = true; }
        }
        if (shut) { r = solve(r.node_pressure, true); ++st.solves; }
    }
    return st;
}

template<class Sys, class Result, class Solve, class Held>
LiftStats allocateLiftGasCore(Sys& system, Result& r, const LiftProblem<typename Sys::ScalarType>& prob,
                              const LiftMethod method, Solve&& solve, Held&& held)
{
    using Scalar = typename Sys::ScalarType;
    LiftStats st;
    const Scalar inc = prob.inc, eco = prob.eco;
    auto oilOf = [&](const Result& x) {
        Scalar oil = 0;
        for (int k = 0; k < system.numWells(); ++k) { oil += system.wells()[k].efficiency * x.well_rate[k]; }
        return oil;
    };
    auto solveFrom = [&](const std::vector<char>& dead) {
        system.restoreDead(dead);
        ++st.solves;
        return solve(r.node_pressure, true);
    };
    // Dead wells: each one's least lift gas that makes it flow (one increment, doubling); the best oil per gas
    // first, again until none pays or the limits stop it -- they share the node and the limits.
    std::vector<char> tried(prob.wells.size(), 0);
    for (bool revived = true; revived;) {
        revived = false;
        std::optional<std::size_t> best_k;
        Scalar best_gain = eco, best_alq = 0;
        Result best_r = r;
        const auto dead = system.committedDead();
        for (std::size_t k = 0; k < prob.wells.size(); ++k) {
            const auto& l = prob.wells[k];
            if (tried[k] || r.well_rate[l.w] > Scalar{0} || held(system.wells()[l.w].name)) { continue; }
            const Scalar a0 = system.wells()[l.w].alq;
            const Scalar from = std::max(a0, l.lo);
            auto start = dead;
            if (static_cast<std::size_t>(l.w) < start.size()) { start[l.w] = 0; }
            auto flows = [&](const Scalar a1) -> std::optional<Result> {
                system.setWellAlq(l.w, a1);
                if (!withinCaps(system, prob)) { return std::nullopt; }
                auto rt = solveFrom(start);
                if (!(rt.converged && rt.well_rate[l.w] > Scalar{0})) { return std::nullopt; }
                return rt;
            };
            // Doubling to the first amount that flows, then halving back to the least on the increments.
            Scalar lo_n = 0, hi_n = 0;      // in increments: lo_n does not flow, hi_n does
            std::optional<Result> at_hi;
            for (Scalar n = 1; from + n * inc <= l.hi * (1 + 1e-9); n *= 2) {
                if ((at_hi = flows(from + n * inc))) { hi_n = n; break; }
                lo_n = n;
            }
            while (at_hi && hi_n - lo_n > 1) {
                const Scalar mid = std::floor((lo_n + hi_n) / 2);
                if (auto rt = flows(from + mid * inc)) { hi_n = mid; at_hi = rt; } else { lo_n = mid; }
            }
            if (at_hi) {
                const Scalar a1 = from + hi_n * inc;
                const Scalar gain = (oilOf(*at_hi) - oilOf(r)) / (a1 - a0);
                if (gain > best_gain) { best_gain = gain; best_k = k; best_alq = a1; best_r = *at_hi; }
            }
            system.setWellAlq(l.w, a0);
            system.restoreDead(dead);
        }
        if (best_k) {
            const int w = prob.wells[*best_k].w;
            system.setWellAlq(w, best_alq);
            auto start = dead;
            if (static_cast<std::size_t>(w) < start.size()) { start[w] = 0; }
            system.restoreDead(start);
            r = solveFrom(start);
            tried[*best_k] = 1;
            ++st.moves;
            revived = r.converged;
        } else {
            (void)system.reducedResidual(r.node_pressure);
        }
    }
    std::vector<int> wells;
    for (const auto& l : prob.wells) { wells.push_back(l.w); }
    if (wells.empty() || !r.converged) { return st; }

    // One round of trials: an increment up and down per well, a route solve each; downs below the minimum,
    // then ups best first within the limits. Repeated until a round moves nothing -- a local optimum on the
    // increments; it is the trials method, and the polish after the gradient methods.
    auto trialRound = [&]() {
        struct Move { std::size_t k; Scalar gradient; };
        std::vector<Move> ups, downs;
        const auto dead = system.committedDead();
        const Scalar base = oilOf(r);
        for (std::size_t k = 0; k < prob.wells.size(); ++k) {
            const auto& l = prob.wells[k];
            if (held(system.wells()[l.w].name)) { continue; }
            const Scalar a0 = system.wells()[l.w].alq;
            for (const int dir : {1, -1}) {
                const Scalar a1 = a0 + dir * inc;
                if (a1 > l.hi * (1 + 1e-9) || a1 < l.lo - Scalar{1e-12}) { continue; }
                system.setWellAlq(l.w, a1);
                auto start = dead;
                if (dir > 0 && static_cast<std::size_t>(l.w) < start.size()) { start[l.w] = 0; }
                const auto rt = solveFrom(start);
                system.setWellAlq(l.w, a0);
                if (!rt.converged) { continue; }
                (dir > 0 ? ups : downs).push_back({k, l.wf * dir * (oilOf(rt) - base) / inc});
            }
        }
        system.restoreDead(dead);
        std::vector<Scalar> before(prob.wells.size());
        for (std::size_t k = 0; k < before.size(); ++k) { before[k] = system.wells()[prob.wells[k].w].alq; }
        int moved = 0;
        for (const auto& m : downs) {
            if (m.gradient < eco) { system.setWellAlq(prob.wells[m.k].w, system.wells()[prob.wells[m.k].w].alq - inc); ++moved; }
        }
        std::sort(ups.begin(), ups.end(), [](const Move& a, const Move& b) { return a.gradient > b.gradient; });
        auto start = dead;
        for (const auto& m : ups) {
            const int w = prob.wells[m.k].w;
            const bool lowered = std::any_of(downs.begin(), downs.end(), [&](const Move& d) { return d.k == m.k && d.gradient < eco; });
            if (m.gradient <= eco || lowered) { continue; }
            system.setWellAlq(w, system.wells()[w].alq + inc);
            if (!withinCaps(system, prob)) { system.setWellAlq(w, system.wells()[w].alq - inc); continue; }
            if (static_cast<std::size_t>(w) < start.size()) { start[w] = 0; }
            ++moved;
        }
        if (moved > 0) {
            const Scalar j0 = liftObjective(system, r, prob);
            auto saved = r;
            r = solveFrom(start);
            // Moves judged one at a time can overshoot together: keep the round only if J rose.
            if (!r.converged || liftObjective(system, r, prob) <= j0) {
                for (std::size_t k = 0; k < before.size(); ++k) { system.setWellAlq(prob.wells[k].w, before[k]); }
                system.restoreDead(dead);
                r = saved;
                (void)system.reducedResidual(r.node_pressure);
                return 0;
            }
        }
        st.moves += moved;
        return moved;
    };
    auto polish = [&]() {
        for (int round = 0; round < 20 && r.converged; ++round) {
            if (trialRound() == 0) { break; }
        }
    };

    if (method == LiftMethod::Trials) {
        polish();
        return st;
    }

    if (method == LiftMethod::Projected) {
        // Gradient ascent on J = oil - eco * lift, projected onto the bounds and limits, with a line search on J;
        // the first trial moves the steepest well four increments.
        std::vector<Scalar> a(prob.wells.size());
        for (std::size_t k = 0; k < a.size(); ++k) { a[k] = system.wells()[prob.wells[k].w].alq; }
        Scalar j = liftObjective(system, r, prob);
        for (int it = 0; it < 20; ++it) {
            const auto g = liftGradients(system, r.node_pressure, wells, Scalar{1e-2} * inc, Scalar{1e-3} * unit::barsa);
            ++st.gradients;
            if (!g.ok) { break; }
            std::vector<Scalar> d(a.size());
            Scalar dmax = 0;
            for (std::size_t k = 0; k < a.size(); ++k) {
                const auto& l = prob.wells[k];
                d[k] = held(system.wells()[l.w].name)
                    ? Scalar{0} : l.wf * g.doil[k] / std::max(Scalar{1} + l.gf * g.dgas[k], Scalar{1e-12}) - eco;
                dmax = std::max(dmax, std::abs(d[k]));
            }
            if (!(dmax > Scalar{0})) { break; }
            const auto dead = system.committedDead();
            bool accepted = false;
            Scalar moved = 0;
            for (Scalar alpha = 4 * inc / dmax; alpha * dmax >= inc / 4; alpha /= 2) {
                auto trial = a;
                for (std::size_t k = 0; k < a.size(); ++k) { trial[k] = a[k] + alpha * d[k]; }
                projectLift(trial, prob);
                auto start = dead;
                for (std::size_t k = 0; k < a.size(); ++k) {
                    system.setWellAlq(prob.wells[k].w, trial[k]);
                    if (trial[k] > a[k] && static_cast<std::size_t>(prob.wells[k].w) < start.size()) { start[prob.wells[k].w] = 0; }
                }
                const auto rt = solveFrom(start);
                if (rt.converged && liftObjective(system, rt, prob) > j + Scalar{1e-9} * std::abs(j)) {
                    for (std::size_t k = 0; k < a.size(); ++k) { moved = std::max(moved, std::abs(trial[k] - a[k])); }
                    a = trial;
                    r = rt;
                    j = liftObjective(system, r, prob);
                    accepted = true;
                    ++st.moves;
                    break;
                }
                for (std::size_t k = 0; k < a.size(); ++k) { system.setWellAlq(prob.wells[k].w, a[k]); }
                system.restoreDead(dead);
            }
            if (!accepted || moved < inc / 2) { break; }
        }
        // To LIFTOPT's increments: down to the lattice, then up where the step still pays.
        for (std::size_t k = 0; k < a.size(); ++k) {
            const auto& l = prob.wells[k];
            const Scalar on = l.lo + inc * std::floor((a[k] - l.lo) / inc + Scalar{1e-6});
            system.setWellAlq(l.w, std::clamp(on, l.lo, l.hi));
        }
        r = solveFrom(system.committedDead());
        polish();
        return st;
    }

    // Steps: each well toward its marginal meeting the minimum, doubling while the direction holds, halving
    // when it flips, down to the increment; the limits after each update, lowest gradient gives back first.
    std::vector<Scalar> step(prob.wells.size(), inc);
    std::vector<int> last(prob.wells.size(), 0);
    for (int inner = 0; inner < 30 && r.converged; ++inner) {
        const auto g = liftGradients(system, r.node_pressure, wells, Scalar{1e-2} * inc, Scalar{1e-3} * unit::barsa);
        ++st.gradients;
        if (!g.ok) { break; }
        bool moved = false;
        for (std::size_t k = 0; k < prob.wells.size(); ++k) {
            const auto& l = prob.wells[k];
            if (last[k] == 2 || held(system.wells()[l.w].name)) { continue; }
            const Scalar alq = system.wells()[l.w].alq;
            const Scalar grad = l.wf * g.doil[k] / std::max(Scalar{1} + l.gf * g.dgas[k], Scalar{1e-12});
            const int dir = grad > eco ? 1 : (grad < eco && alq > l.lo) ? -1 : 0;
            if (dir == 0) { continue; }
            if (last[k] != 0 && dir != last[k]) {
                step[k] = std::max(inc, inc * std::floor(step[k] / (2 * inc)));
                if (step[k] <= inc) {
                    if (dir < 0) { system.setWellAlq(l.w, std::max(l.lo, alq - inc)); moved = true; }
                    last[k] = 2;
                    continue;
                }
            } else if (last[k] == dir) {
                step[k] = std::min(2 * step[k], std::max(l.hi - l.lo, inc));
            }
            const Scalar target = std::clamp(alq + dir * step[k], l.lo, l.hi);
            if (target == alq) { continue; }
            system.setWellAlq(l.w, target);
            last[k] = dir;
            moved = true;
        }
        for (const auto& c : prob.caps) {
            std::vector<std::pair<Scalar, std::size_t>> under;
            Scalar total = 0;
            for (std::size_t k = 0; k < prob.wells.size(); ++k) {
                if (std::find(c.wells.begin(), c.wells.end(), prob.wells[k].w) == c.wells.end()) { continue; }
                total += system.wells()[prob.wells[k].w].alq;
                under.emplace_back(prob.wells[k].wf * g.doil[k], k);
            }
            std::sort(under.begin(), under.end());
            for (const auto& [grad, k] : under) {
                if (total <= c.cap * (1 + 1e-9)) { break; }
                const int w = prob.wells[k].w;
                const Scalar cut = std::min(system.wells()[w].alq - prob.wells[k].lo, inc * std::ceil((total - c.cap) / inc - 1e-9));
                if (cut <= 0) { continue; }
                system.setWellAlq(w, system.wells()[w].alq - cut);
                total -= cut;
                moved = true;
            }
        }
        if (!moved) { break; }
        ++st.moves;
        auto start = system.committedDead();
        for (const auto& l : prob.wells) {
            if (system.wells()[l.w].alq > 0 && static_cast<std::size_t>(l.w) < start.size()) { start[l.w] = 0; }
        }
        r = solveFrom(start);
    }
    polish();
    return st;
}

/// The yardstick: every allocation on LIFTOPT's increments within the bounds and limits, a route solve each,
/// the best J. Wells start alive; at most max_points allocations.
template<class Sys, class Result, class Solve>
std::pair<typename Sys::ScalarType, std::vector<typename Sys::ScalarType>>
bestLiftByEnumeration(Sys& system, const Result& r, const LiftProblem<typename Sys::ScalarType>& prob, Solve&& solve,
                      const long max_points = 200000)
{
    using Scalar = typename Sys::ScalarType;
    const auto dead0 = system.committedDead();
    auto alive = dead0;
    for (const auto& l : prob.wells) { if (static_cast<std::size_t>(l.w) < alive.size()) { alive[l.w] = 0; } }
    std::vector<Scalar> a(prob.wells.size()), best_a;
    Scalar best = -std::numeric_limits<Scalar>::max();
    long points = 0;
    std::function<void(std::size_t)> walk = [&](const std::size_t k) {
        if (points >= max_points) { return; }
        if (k == prob.wells.size()) {
            if (!withinCaps(system, prob)) { return; }
            ++points;
            system.restoreDead(alive);
            const auto rt = solve(r.node_pressure, true);
            if (rt.converged) {
                const Scalar j = liftObjective(system, rt, prob);
                if (j > best) { best = j; best_a = a; }
            }
            return;
        }
        const auto& l = prob.wells[k];
        for (Scalar v = l.lo; v <= l.hi * (1 + 1e-9); v += prob.inc) {
            a[k] = v;
            system.setWellAlq(l.w, v);
            if (!withinCaps(system, prob)) { break; }
            walk(k + 1);
        }
        system.setWellAlq(l.w, l.lo);
    };
    walk(0);
    system.restoreDead(dead0);
    return {best, best_a};
}

/// The lift problem as the controller dumps it next to the system: "gaslift inc eco", "lifted w lo hi wf gf",
/// "liftcap cap w...".
template<class Scalar>
void writeLiftProblem(const LiftProblem<Scalar>& prob, std::ostream& os)
{
    os.precision(17);
    os << "gaslift " << prob.inc << ' ' << prob.eco << '\n';
    for (const auto& l : prob.wells) { os << "lifted " << l.w << ' ' << l.lo << ' ' << l.hi << ' ' << l.wf << ' ' << l.gf << '\n'; }
    for (const auto& c : prob.caps) {
        os << "liftcap " << c.cap;
        for (const int w : c.wells) { os << ' ' << w; }
        os << '\n';
    }
}

template<class Scalar>
LiftProblem<Scalar> readLiftProblem(std::istream& is)
{
    LiftProblem<Scalar> prob;
    std::string line;
    while (std::getline(is, line)) {
        std::istringstream in(line);
        std::string tag;
        if (!(in >> tag)) { continue; }
        if (tag == "gaslift") { in >> prob.inc >> prob.eco; }
        else if (tag == "lifted") { LiftWell<Scalar> l{}; in >> l.w >> l.lo >> l.hi >> l.wf >> l.gf; prob.wells.push_back(l); }
        else if (tag == "liftcap") { LiftCap<Scalar> c{}; in >> c.cap; int w; while (in >> w) { c.wells.push_back(w); } prob.caps.push_back(c); }
    }
    return prob;
}

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_GAS_LIFT_HEADER_INCLUDED
