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
#ifndef OPM_NETWORK_REDUCED_SOLVE_HEADER_INCLUDED
#define OPM_NETWORK_REDUCED_SOLVE_HEADER_INCLUDED

#include <opm/simulators/wells/network/NetworkSolve.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

namespace Opm::NetworkSolve {

template<class Scalar>
struct ReducedResult
{
    bool converged = false;
    int iterations = 0;
    int evaluations = 0;      // residual evaluations, differences included
    int stalls = 0;           // steps the line search could not improve on
    int cliffs = 0;           // steps cut back to keep a well alive
    /// Settled on the alive side of a cliff the residual could not resolve:
    /// a well that dies when the pressure rises and revives when it falls
    /// has no fixed point, and the answer taken is the one with it flowing.
    bool on_cliff = false;
    int held_at_cliff = 0;    // wells given their cliff rate instead of dying
    int set_changes = 0;      // iterations after which the tree walk chose differently
    Scalar residual = 0;
    int off_axis = 0;         // lookups the answer needed off a table axis
    std::vector<Scalar> node_pressure;
    std::vector<Scalar> well_rate;
    std::string sets;
};

/// The reduced Jacobian by elimination instead of differences. At the state
/// the reduced residual was evaluated at, the full system's rows for the
/// same active set are all satisfied except the node rows, and their
/// assembled Jacobian J splits into node pressures p and the rest x; the
/// derivative of the node rows along the manifold the rest defines is the
/// Schur complement J_pp - J_px J_xx^-1 J_xp. One assembly and N solves of
/// the well-and-group block, no table lookups beyond the gradients the
/// assembly already takes.
template<class Sys>
DenseMatrix<typename Sys::ScalarType>
reducedJacobianByElimination(const Sys& system)
{
    using Scalar = typename Sys::ScalarType;
    const auto J = system.jacobian(system.reducedState());
    const int n = system.numNodes(), nf = system.size(), nx = nf - n;
    DenseMatrix<Scalar> Jxx(nx);
    for (int i = 0; i < nx; ++i) {
        for (int k = 0; k < nx; ++k) { Jxx(i, k) = J(n + i, n + k); }
    }
    DenseMatrix<Scalar> S(n);
    std::vector<Scalar> col(nx), y;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < nx; ++i) { col[i] = J(n + i, j); }
        if (!Jxx.solve(col, y)) { y.assign(nx, Scalar{0}); }
        for (int i = 0; i < n; ++i) {
            Scalar s = J(i, j);
            for (int k = 0; k < nx; ++k) { s -= J(i, n + k) * y[k]; }
            S(i, j) = s;
        }
    }
    return S;
}

/// The same step without forming the Schur complement: with every row but
/// the node rows satisfied at the reduced state, J dz = -[r; 0] has the
/// reduced Newton step as its node-pressure part. One factorisation instead
/// of N.
template<class Sys>
std::vector<typename Sys::ScalarType>
reducedStepByElimination(const Sys& system, const std::vector<typename Sys::ScalarType>& r)
{
    using Scalar = typename Sys::ScalarType;
    const auto J = system.jacobian(system.reducedState());
    const int n = system.numNodes(), nf = system.size();
    std::vector<Scalar> rhs(nf, Scalar{0}), dz;
    for (int i = 0; i < n; ++i) { rhs[i] = -r[i]; }
    if (!J.solve(rhs, dz)) { return {}; }
    return std::vector<Scalar>(dz.begin(), dz.begin() + n);
}

/// Newton on the node pressures alone, with the tree walk as the well model:
/// r(p) = p - branch(p_up, q(c(p))), c the wells' capacities at p, q what the
/// walk makes of them. Continuous and piecewise smooth; the Jacobian is that
/// of the piece the iterate is on -- by differences, or by elimination from
/// the full system's -- and a backtracking line search takes the steps that
/// cross into another.
/// What a well does at a cliff the residual asks to cross: die (the well
/// model's rule, sticky for the solve) or hold the rate it had just before
/// (a wellhead choke).
enum class CliffRule { Die, Hold };

template<class Sys>
ReducedResult<typename Sys::ScalarType>
solveReduced(Sys& system,
             const std::vector<typename Sys::ScalarType>& node_pressure_guess,
             const Parameters<typename Sys::ScalarType> params,
             const bool eliminate = false,
             const CliffRule cliff_rule = CliffRule::Die)
{
    using Scalar = typename Sys::ScalarType;
    ReducedResult<Scalar> out;
    system.flatTargetAsTree();
    system.setGroupActiveSet(true);
    system.setTreeFrozen(false);
    system.setExactPotential(true);
    system.setDeadWhenCannotLift(true);
    system.resetDead();
    system.resetCliffRates();
    const int n = system.numNodes();
    auto p = node_pressure_guess;
    auto norm = [](const std::vector<Scalar>& r) {
        Scalar worst = 0;
        for (const auto e : r) { worst = std::max(worst, std::abs(e)); }
        return worst;
    };
    auto r = system.reducedResidual(p);
    ++out.evaluations;
    // A well that cannot lift at the starting pressures and that the well
    // model already has at zero stays dead; one the well model has flowing
    // is given its chance. Without this a well revives on the first step
    // down and dies on the next one up, forever.
    system.commitDeadNotFlowing();
    bool at_cliff = false;
    static const bool trace = std::getenv("OPM_REDUCED_TRACE") != nullptr;
    auto alive_p = p;
    std::string alive_set;
    std::string last_set = system.treeSignature();
    out.sets = last_set;
    const Scalar max_step = Scalar{50} * unit::barsa, floor = unit::barsa;
    std::vector<std::pair<std::string, Scalar>> recent;   // set and residual, last few iterates
    for (int it = 1; it <= params.max_iterations; ++it) {
        out.iterations = it;
        out.residual = norm(r);
        if (trace) {
            std::string ps;
            for (int i = 1; i <= n; ++i) { ps += fmt::format(" {:.3f}", p[i] / unit::barsa); }
            std::fprintf(stderr, "[reduced] it %d residual %.4g set %s p%s\n", it, out.residual,
                         system.treeSignature().c_str(), ps.c_str());
        }
        if (out.residual < params.tolerance) {
            out.converged = true;
            break;
        }
        // Only two sets seen over the last eight iterates, differing in which
        // wells can lift, every residual nearly converged: a cliff with no
        // fixed point. Take the side with the wells flowing and stop.
        const auto sig = system.treeSignature();
        const auto dead = [](const std::string& x) { return std::count(x.begin(), x.end(), 'S'); };
        recent.emplace_back(sig, out.residual);
        if (recent.size() > 8) { recent.erase(recent.begin()); }
        if (alive_set.empty() || dead(sig) <= dead(alive_set)) { alive_set = sig; alive_p = p; }
        if (recent.size() == 8) {
            std::string a = recent[0].first, b;
            bool two = true;
            Scalar worst = 0;
            for (const auto& [set, res] : recent) {
                worst = std::max(worst, res);
                if (set == a) { continue; }
                if (b.empty()) { b = set; } else if (set != b) { two = false; }
            }
            if (two && !b.empty() && dead(a) != dead(b) && worst < Scalar{50} * params.tolerance) {
                p = alive_p;
                r = system.reducedResidual(p); ++out.evaluations;
                out.residual = norm(r);
                out.converged = true;
                out.on_cliff = true;
                break;
            }
        }
        std::vector<Scalar> dx;
        if (eliminate) {
            dx = reducedStepByElimination(system, r);
            if (dx.empty()) { break; }
        } else {
            DenseMatrix<Scalar> J(n);
            const Scalar h = Scalar{1e-3} * unit::barsa;
            for (int j = 0; j < n; ++j) {
                auto pj = p;
                pj[j + 1] += h;
                const auto rj = system.reducedResidual(pj);
                ++out.evaluations;
                for (int i = 0; i < n; ++i) { J(i, j) = (rj[i] - r[i]) / h; }
            }
            std::vector<Scalar> negative(n);
            for (int i = 0; i < n; ++i) { negative[i] = -r[i]; }
            if (!J.solve(negative, dx)) { break; }
        }
        Scalar alpha = 1;
        for (int i = 0; i < n; ++i) {
            if (std::abs(dx[i]) > max_step) { alpha = std::min(alpha, max_step / std::abs(dx[i])); }
            if (p[i + 1] + alpha * dx[i] < floor) { alpha = std::min(alpha, (floor - p[i + 1]) / dx[i]); }
        }
        auto trial = [&](const Scalar step) {
            auto q = p;
            for (int i = 0; i < n; ++i) { q[i + 1] += step * dx[i]; }
            return q;
        };
        // A cliff -- a well the full step would leave unable to lift -- is
        // crossed only when the residual asks for it twice: the first step
        // that would cross is cut back to the alive side, and if the next
        // step from there still points across, the death is committed and
        // the step taken. A well that dies at a trial point the line search
        // then discards must not die for good.
        const auto dead_here = system.deadNow();
        auto newDeath = [&](const std::vector<char>& d) {
            for (std::size_t w = 0; w < d.size(); ++w) { if (d[w] && !dead_here[w]) { return true; } }
            return false;
        };
        {
            const auto pt = trial(alpha);
            (void)system.reducedResidual(pt);
            ++out.evaluations;
            if (newDeath(system.deadNow())) {
                if (!at_cliff) {
                    // Cut back to the alive side: the largest step with no new death.
                    Scalar lo = 0, hi = alpha;
                    for (int k = 0; k < 10; ++k) {
                        const Scalar mid = Scalar{0.5} * (lo + hi);
                        (void)system.reducedResidual(trial(mid));
                        ++out.evaluations;
                        (newDeath(system.deadNow()) ? hi : lo) = mid;
                    }
                    alpha = lo;
                    at_cliff = true;
                    ++out.cliffs;
                } else if (cliff_rule == CliffRule::Hold) {
                    // Asked twice: the wells that would die hold the rate they
                    // have on the alive side, and the step is taken.
                    (void)system.reducedResidual(p);
                    ++out.evaluations;
                    const auto dying = [&] {
                        (void)system.reducedResidual(pt);
                        ++out.evaluations;
                        return system.deadNow();
                    }();
                    (void)system.reducedResidual(p);
                    ++out.evaluations;
                    for (std::size_t w = 0; w < dying.size(); ++w) {
                        if (dying[w] && !dead_here[w]) {
                            system.setCliffRate(static_cast<int>(w), system.ownAllowance(static_cast<int>(w)));
                            ++out.held_at_cliff;
                        }
                    }
                    at_cliff = false;
                } else {
                    // Asked twice: cross, and the death holds.
                    (void)system.reducedResidual(pt);
                    system.commitDead();
                    at_cliff = false;
                }
            } else {
                at_cliff = false;
            }
        }
        bool accepted = false;
        Scalar step = alpha;
        for (int k = 0; k < 12 && step > Scalar{0}; ++k) {
            const auto pt = trial(step);
            const auto rt = system.reducedResidual(pt);
            ++out.evaluations;
            if (norm(rt) < out.residual || (at_cliff && k == 0)) {
                p = pt; r = rt; accepted = true;
                break;
            }
            step *= Scalar{0.5};
        }
        if (!accepted) {
            // Nothing along the direction improves: a kink between here and
            // there. Take the step anyway and let the next piece's Jacobian
            // say where to go.
            ++out.stalls;
            p = trial(alpha);
            r = system.reducedResidual(p);
            ++out.evaluations;
        }
        const auto set = system.treeSignature();
        if (set != last_set) { ++out.set_changes; out.sets += " " + set; }
        last_set = set;
    }
    out.node_pressure = p;
    out.well_rate = system.wellRates(system.reducedState());
    system.resetOffAxis();
    (void)system.reducedResidual(p);
    out.off_axis = system.offAxisLookups();
    return out;
}

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_REDUCED_SOLVE_HEADER_INCLUDED
