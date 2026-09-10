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
#ifndef OPM_NETWORK_LEGACY_SOLVE_HEADER_INCLUDED
#define OPM_NETWORK_LEGACY_SOLVE_HEADER_INCLUDED

#include <opm/simulators/wells/network/NetworkSolve.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace Opm::NetworkSolve {

/// The simulator's own rules for a production network, as a standalone
/// route on the same system the other routes solve, so the four can be put
/// side by side. Not the simulator's code -- that is bound to the well and
/// group state -- but its rules, read off BlackoilWellModelNetwork_impl.hpp
/// and BlackoilWellModelNetworkGeneric.cpp, with the simulator's defaults:
///
///   - every well at its dynamic thp limit, the node pressure, takes the
///     most restrictive of its tubing crossing, its bhp limit and its rate
///     limit; a well the well model found dead above some pressure is dead;
///   - an autochoke group's thp is found by brute force: a 1 bar sweep up
///     from the terminal pressure for a sign change of (rate - target),
///     300 samples across that bracket, then bisection; if nothing brackets,
///     the choke is open; the group thp is never below the node above;
///   - node pressures follow from the rates through the branch tables and
///     are moved a damped, capped step toward them (0.1, 5 bar), up to 100
///     times, until the largest move is under the balance tolerance.
///
/// Guide-rate shares under a plain group target are not here; none of the
/// cases this runs on has one that binds without a choke.
template<class Scalar>
struct LegacyParameters
{
    Scalar damping = 0.1;
    Scalar max_update = 5.0 * unit::barsa;
    Scalar tolerance = 0.01 * unit::barsa;       // NETBALAN item 2
    Scalar thp_tolerance = 0.01 * unit::barsa;   // NETBALAN item 4
    int max_iterations = 100;                    // NetworkMaxSubIterations
    int bracket_samples = 300;                   // NetworkAutochokeBracketSamples
};

template<class Scalar>
struct LegacyResult
{
    bool converged = false;
    int iterations = 0;
    Scalar imbalance = 0;
    std::vector<Scalar> node_pressure;
    std::vector<Scalar> well_rate;
    std::string controls;               // one letter per well: T, B, O, D (dead), S
    std::optional<Scalar> choke_thp;    // the group thp found, if a choke is held
    bool choke_bracketed = false;
    Scalar choke_mismatch = 0;          // (rate - target)/target at the thp taken
    int choke_evaluations = 0;
    int off_axis = 0;                   // lookups the answer needed off a table axis
};

template<class Sys>
LegacyResult<typename Sys::ScalarType>
solveLegacy(Sys& system,
            const std::vector<typename Sys::ScalarType>& node_pressure_guess,
            const LegacyParameters<typename Sys::ScalarType> lp = {})
{
    using Scalar = typename Sys::ScalarType;
    constexpr int NP = Sys::NP;
    LegacyResult<Scalar> out;
    system.setExactPotential(true);
    const int N = system.numNodes(), W = system.numWells();
    auto p = node_pressure_guess;
    p.resize(N + 1);
    p[0] = system.terminalPressure();

    // A well solved at a dynamic thp limit, with switching.
    auto rateAt = [&](const int w, const Scalar thp, char& letter) {
        const auto& well = system.wells()[w];
        if (well.shut) { letter = 'S'; return Scalar{0}; }
        if (well.pinned) { letter = 'O'; return well.oil_rate_limit; }
        Scalar q = std::max(system.ipr(well, 1, well.bhp_limit), Scalar{0});
        letter = 'B';
        if (well.oil_rate_limit > Scalar{0} && well.oil_rate_limit < q) { q = well.oil_rate_limit; letter = 'O'; }
        if (well.vfp_table > 0) {
            if (well.dead_above > Scalar{0} && thp >= well.dead_above) { letter = 'D'; return Scalar{0}; }
            const Scalar c = system.thpPotential(well, thp);
            if (!(c > Scalar{0})) { letter = 'D'; return Scalar{0}; }
            if (c < q) { q = c; letter = 'T'; }
        }
        return q;
    };
    auto phases = [&](const int w, const Scalar q_oil) {
        const auto& well = system.wells()[w];
        std::array<Scalar, NP> q{};
        if (!(q_oil > Scalar{0}) || !(well.ipr_b[1] < Scalar{0})) { return q; }
        const Scalar bhp = (q_oil - well.ipr_a[1]) / well.ipr_b[1];
        for (int ph = 0; ph < NP; ++ph) { q[ph] = std::max(system.ipr(well, ph, bhp), Scalar{0}); }
        return q;
    };

    // Chokes: one per node with a target; the wells at it are the group's.
    int choke = -1;
    for (int n = 1; n <= N; ++n) { if (system.isChoke(n)) { choke = n; } }
    auto upstreamOf = [&](const int n) { return system.nodes()[n].parent; };

    std::vector<Scalar> q_oil(W, Scalar{0});
    std::string letters(W, '-');
    std::optional<Scalar> autochoke_thp;
    for (int it = 1; it <= lp.max_iterations; ++it) {
        out.iterations = it;
        // The choke's thp, from the wells beneath it.
        Scalar group_thp = 0;
        if (choke >= 0) {
            const Scalar target = system.chokeTarget(choke);
            const Scalar nodal = p[upstreamOf(choke)];
            auto mismatch = [&](const Scalar thp) {
                ++out.choke_evaluations;
                Scalar sum = 0;
                for (int w = 0; w < W; ++w) {
                    if (system.wells()[w].node != choke) { continue; }
                    char l; sum += system.wells()[w].efficiency * rateAt(w, thp, l);
                }
                return (sum - target) / target;
            };
            std::array<Scalar, 2> range{};
            bool have_range = false;
            if (!autochoke_thp) {
                // The 1 bar sweep up from the terminal pressure.
                Scalar lo = system.terminalPressure(), hi = lo;
                Scalar f_lo = mismatch(lo);
                for (int i = 1; i <= 1000; ++i) {
                    hi = system.terminalPressure() + Scalar{1e5} * i;
                    const Scalar f_hi = mismatch(hi);
                    if (f_hi * f_lo <= Scalar{0}) { lo = hi - Scalar{1e5}; have_range = true; break; }
                    f_lo = f_hi;
                }
                if (have_range) { range = {Scalar{0.9} * lo, Scalar{1.1} * hi}; }
            } else {
                range = {Scalar{0.9} * *autochoke_thp, Scalar{1.1} * *autochoke_thp};
                have_range = true;
            }
            std::optional<Scalar> found;
            if (have_range && (!autochoke_thp || *autochoke_thp > nodal)) {
                // The fine bracket, then bisection to the thp tolerance.
                const Scalar step = (range[1] - range[0]) / lp.bracket_samples;
                Scalar lo = range[0], f_lo = mismatch(lo), hi = lo, f_hi = f_lo;
                bool bracketed = false;
                for (int i = 0; i <= lp.bracket_samples; ++i) {
                    hi = range[0] + step * i;
                    f_hi = mismatch(hi);
                    if (f_hi * f_lo <= Scalar{0}) { bracketed = true; break; }
                    lo = hi; f_lo = f_hi;
                }
                out.choke_bracketed = bracketed;
                if (bracketed) {
                    for (int k = 0; k < 100 && (hi - lo) > lp.thp_tolerance; ++k) {
                        const Scalar mid = Scalar{0.5} * (lo + hi);
                        const Scalar f_mid = mismatch(mid);
                        if (f_mid * f_lo <= Scalar{0}) { hi = mid; f_hi = f_mid; } else { lo = mid; f_lo = f_mid; }
                    }
                    found = Scalar{0.5} * (lo + hi);
                }
            } else if (autochoke_thp) {
                found = autochoke_thp;
            }
            if (found) { autochoke_thp = found; }
            group_thp = autochoke_thp ? std::max(*autochoke_thp, nodal) : nodal;
            out.choke_thp = autochoke_thp;
            out.choke_mismatch = mismatch(group_thp);
        }
        // Every well at its node pressure -- the choke's wells at the group thp.
        for (int w = 0; w < W; ++w) {
            const int n = system.wells()[w].node;
            const Scalar thp = (n == choke && choke >= 0) ? group_thp : p[n];
            q_oil[w] = rateAt(w, thp, letters[w]);
        }
        // Pressures from the rates: sums up the tree, tables down it.
        std::vector<std::array<Scalar, NP>> Q(N + 1);
        for (int n = N; n >= 1; --n) {
            auto q = system.nodeSource(n);
            for (int w = 0; w < W; ++w) {
                if (system.wells()[w].node != n) { continue; }
                const auto qw = phases(w, q_oil[w]);
                for (int ph = 0; ph < NP; ++ph) { q[ph] += system.wells()[w].efficiency * qw[ph]; }
                q[2] += system.wells()[w].efficiency * system.wells()[w].lift_gas;
            }
            for (int c = 1; c <= N; ++c) {
                if (system.nodes()[c].parent != n) { continue; }
                for (int ph = 0; ph < NP; ++ph) { q[ph] += system.nodes()[c].efficiency * Q[c][ph]; }
            }
            Q[n] = q;
        }
        std::vector<Scalar> computed(N + 1, system.terminalPressure());
        for (int n = 1; n <= N; ++n) {
            const auto& node = system.nodes()[n];
            const Scalar up = computed[node.parent];
            computed[n] = (n == choke) ? std::max(group_thp, up)
                : (node.vfp_table != NoTable) ? system.tableBhp(node.vfp_table, up, Q[n], system.branchAlq(n)) : up;
        }
        // The damped, capped move.
        Scalar imbalance = 0;
        for (int n = 1; n <= N; ++n) {
            const Scalar d = computed[n] - p[n];
            imbalance = std::max(imbalance, std::abs(d));
            const Scalar move = std::clamp(lp.damping * d, -lp.max_update, lp.max_update);
            p[n] += move;
        }
        out.imbalance = imbalance;
        if (imbalance < lp.tolerance) {
            out.converged = true;
            // The same count the other routes report: the answer's lookups.
            system.resetOffAxis();
            const typename Sys::CountScope counting(system, true);
            for (int w = 0; w < W; ++w) {
                const int n = system.wells()[w].node;
                char l;
                (void)rateAt(w, (n == choke && choke >= 0) ? group_thp : p[n], l);
            }
            for (int n = 1; n <= N; ++n) {
                const auto& node = system.nodes()[n];
                if (node.vfp_table != NoTable && n != choke) {
                    (void)system.tableBhp(node.vfp_table, p[node.parent], Q[n], system.branchAlq(n));
                }
            }
            out.off_axis = system.offAxisLookups();
            break;
        }
    }
    out.node_pressure = p;
    out.well_rate = q_oil;
    out.controls = letters;
    return out;
}

/// The legacy update with the walk as its well model: the reduced residual
/// gives the pressures the current rates imply, and each node moves a
/// damped, capped step toward them -- the simulator's fixed point, on a
/// system with a group tree.
template<class Sys>
LegacyResult<typename Sys::ScalarType>
solveLegacyOnWalk(Sys& system,
                  const std::vector<typename Sys::ScalarType>& node_pressure_guess,
                  const LegacyParameters<typename Sys::ScalarType> lp = {})
{
    using Scalar = typename Sys::ScalarType;
    LegacyResult<Scalar> out;
    system.flatTargetAsTree();
    system.setGroupActiveSet(true);
    system.setTreeFrozen(false);
    system.setExactPotential(true);
    system.setDeadWhenCannotLift(true);
    system.resetDead();
    system.resetCliffRates();
    const int N = system.numNodes();
    auto p = node_pressure_guess;
    p.resize(N + 1);
    p[0] = system.terminalPressure();
    for (int it = 1; it <= lp.max_iterations; ++it) {
        out.iterations = it;
        const auto r = system.reducedResidual(p);       // r_n = (p_n - computed_n) / scale
        Scalar imbalance = 0;
        for (int n = 1; n <= N; ++n) {
            const Scalar d = -r[n - 1] * unit::barsa;   // computed - applied
            imbalance = std::max(imbalance, std::abs(d));
            p[n] += std::clamp(lp.damping * d, -lp.max_update, lp.max_update);
        }
        out.imbalance = imbalance;
        if (imbalance < lp.tolerance) { out.converged = true; break; }
    }
    (void)system.reducedResidual(p);
    out.node_pressure = p;
    out.well_rate = system.wellRates(system.reducedState());
    out.controls.clear();
    for (int w = 0; w < system.numWells(); ++w) { out.controls += system.controlLetter(w); }
    return out;
}

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_LEGACY_SOLVE_HEADER_INCLUDED
