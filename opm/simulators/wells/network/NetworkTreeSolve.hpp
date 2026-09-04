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
#ifndef OPM_NETWORK_TREE_SOLVE_HEADER_INCLUDED
#define OPM_NETWORK_TREE_SOLVE_HEADER_INCLUDED

#include <opm/simulators/wells/network/NetworkSolve.hpp>

#include <algorithm>
#include <string>
#include <vector>

namespace Opm::NetworkSolve {

/// What to do when the tree pass hands back a set that was solved already.
enum class OnRevisit { Stop, FischerBurmeister };

template<class Scalar>
struct TreeResult
{
    Result<Scalar> result;            // the last inner solve
    int passes = 0;                   // tree resolutions
    int inner_iterations = 0;         // Newton iterations over all passes
    /// The set the tree pass picks at the answer is the set it was solved with.
    bool consistent = false;
    bool cycled = false;
    bool fell_back = false;
    std::vector<std::string> sets;    // one signature per pass
};

/// Network and group tree together, the tree's active set decided outside
/// the Newton.
///
/// Each pass measures every well's capacity at the current node pressures --
/// the same IPR and tubing curve the Newton uses -- walks the tree with them
/// (resolveTree, which is Stein's balancer used for the set alone), freezes
/// the result and solves the network with it. The pass after a converged
/// solve re-walks the tree at the converged pressures: if nothing moves, the
/// balancer's allocation at those pressures is a solution of the network,
/// and so of the whole system. That is the stopping test, and it is the only
/// one. A set that comes back a second time is a cycle, and is not solved
/// again.
template<class Sys, class Globalisation>
TreeResult<typename Sys::ScalarType>
solveWithTree(Sys& system,
              const std::vector<typename Sys::ScalarType>& node_pressure_guess,
              const Parameters<typename Sys::ScalarType> params,
              Globalisation globalisation,
              const int max_passes = 10,
              const OnRevisit on_revisit = OnRevisit::Stop)
{
    using Scalar = typename Sys::ScalarType;
    TreeResult<Scalar> out;
    system.setGroupActiveSet(true);
    system.setTreeFrozen(true);
    auto p = node_pressure_guess;
    auto x = system.start(p);
    for (int pass = 1; pass <= max_passes; ++pass) {
        out.passes = pass;
        system.updateControls(x);
        const bool moved = system.resolveTree(x);
        const auto set = system.treeSignature();
        if (pass > 1 && !moved && out.result.converged) {
            out.consistent = true;
            break;
        }
        if (std::find(out.sets.begin(), out.sets.end(), set) != out.sets.end()) {
            out.cycled = true;
            break;
        }
        out.sets.push_back(set);
        out.result = solve(system, p, params, globalisation);
        out.inner_iterations += out.result.iterations;
        if (!out.result.converged) {
            break;
        }
        p = out.result.node_pressure;
        x = out.result.state;
    }
    if (out.cycled && on_revisit == OnRevisit::FischerBurmeister) {
        // The continuous formulation, from where the cycle left off.
        system.setGroupActiveSet(false);
        system.setTreeFrozen(false);
        out.result = solve(system, p, params, globalisation);
        out.inner_iterations += out.result.iterations;
        out.fell_back = true;
    }
    system.setTreeFrozen(false);
    return out;
}

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_TREE_SOLVE_HEADER_INCLUDED
