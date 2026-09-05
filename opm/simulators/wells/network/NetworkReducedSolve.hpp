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

#include <algorithm>
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
    int set_changes = 0;      // iterations after which the tree walk chose differently
    Scalar residual = 0;
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
template<class Sys>
ReducedResult<typename Sys::ScalarType>
solveReduced(Sys& system,
             const std::vector<typename Sys::ScalarType>& node_pressure_guess,
             const Parameters<typename Sys::ScalarType> params,
             const bool eliminate = false)
{
    using Scalar = typename Sys::ScalarType;
    ReducedResult<Scalar> out;
    system.setGroupActiveSet(true);
    system.setTreeFrozen(false);
    system.setExactPotential(true);
    const int n = system.numNodes();
    auto p = node_pressure_guess;
    auto norm = [](const std::vector<Scalar>& r) {
        Scalar worst = 0;
        for (const auto e : r) { worst = std::max(worst, std::abs(e)); }
        return worst;
    };
    auto r = system.reducedResidual(p);
    ++out.evaluations;
    std::string last_set = system.treeSignature();
    out.sets = last_set;
    const Scalar max_step = Scalar{50} * unit::barsa, floor = unit::barsa;
    for (int it = 1; it <= params.max_iterations; ++it) {
        out.iterations = it;
        out.residual = norm(r);
        if (out.residual < params.tolerance) {
            out.converged = true;
            break;
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
        bool accepted = false;
        Scalar step = alpha;
        for (int k = 0; k < 12; ++k) {
            const auto pt = trial(step);
            const auto rt = system.reducedResidual(pt);
            ++out.evaluations;
            if (norm(rt) < out.residual) {
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
    return out;
}

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_REDUCED_SOLVE_HEADER_INCLUDED
