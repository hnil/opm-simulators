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
#ifndef OPM_NETWORK_JUDGE_HEADER_INCLUDED
#define OPM_NETWORK_JUDGE_HEADER_INCLUDED

#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/simulators/wells/network/NetworkProductionSystem.hpp>

#include <fmt/format.h>

#include <array>
#include <cmath>
#include <map>
#include <string>
#include <vector>

// The judge: scores an answer -- node pressures, oil rates, the control of
// each well -- against every row of a production system, so that any route's
// answer can be checked independently of how it was found.

namespace Opm::NetworkSolve {

using Opm::unit::barsa;

/// hysteresis: wells shut although a crossing exists at the settled
/// pressure -- the cliff's answer, shut stays shut; not a violation of
/// the rows, but not a strict solution either, so counted apart.
struct Verdict
{
    bool ok = true;
    int hysteresis = 0;
    std::map<std::string, int> violations;
    std::vector<int> groups_over;        // groups above their target with a well that could be reduced
    std::vector<double> over_ratio;      // their rate on mode over the target
    int over_irreducible = 0;            // groups above their target with nothing left to reduce
    std::vector<int> wells_unliftable;   // non-thp wells whose tubing cannot lift their rate
};

template<class Sys>
Verdict verifyAnswer(const Sys& sys, const std::vector<double>& p, const std::vector<double>& q_oil,
                     const std::string& controls, const double dp_tol = 0.05 * barsa, const double r_tol = 0.005)
{
    constexpr int NP = Sys::NP;
    Verdict v;
    auto fail = [&](const std::string& what) { v.ok = false; ++v.violations[what]; };
    std::vector<std::array<double, NP>> q(sys.numWells());
    std::vector<double> bhp(sys.numWells());
    double in_group = 0.0;
    bool any_held = false;
    for (int w = 0; w < sys.numWells(); ++w) {
        const auto& well = sys.wells()[w];
        const char c = controls[w];
        const double qo = q_oil[w];
        // The pressure at the well head: its own thp limit on the no-network route.
        const double pw = well.own_thp > 0.0 ? well.own_thp : p[well.node];
        bhp[w] = (qo > 0.0 && well.ipr_b[1] < 0.0) ? (qo - well.ipr_a[1]) / well.ipr_b[1] : well.bhp_limit;
        for (int ph = 0; ph < NP; ++ph) { q[w][ph] = qo > 0.0 ? std::max(sys.ipr(well, ph, bhp[w]), 0.0) : 0.0; }
        if (well.shut || c == 'S' || c == 'D') {
            if (qo > 1e-9) { fail("shut well producing"); }
            // The branch rule: shut only because nothing can lift it.
            if (!well.shut && well.vfp_table > 0 && sys.thpPotential(well, pw) > 0.0
                && !(well.dead_above > 0.0 && pw >= well.dead_above)) { ++v.hysteresis; }
            continue;
        }
        if (!(qo > 0.0)) { fail("open well at zero"); continue; }
        if (well.pinned) {
            if (std::abs(qo - well.oil_rate_limit) > r_tol * well.oil_rate_limit) { fail("pinned well off its rate"); }
            continue;
        }
        if (well.oil_rate_limit > 0.0 && qo > well.oil_rate_limit * (1 + r_tol)) { fail("above own rate limit"); }
        if (bhp[w] < well.bhp_limit - dp_tol) { fail("below bhp limit"); }
        if (well.vfp_table > 0) {
            const double need = sys.tableBhp(well, pw, q[w]) - well.vfp_dp;
            if (c == 'T') {
                if (std::abs(bhp[w] - need) > dp_tol) { fail("thp well off its tubing curve"); }
                if (sys.thpPotential(well, pw) > 0.0
                    && std::abs(qo - sys.thpPotential(well, pw)) > r_tol * qo) { fail("thp well not on the stable crossing"); }
            } else if (need > bhp[w] + dp_tol && qo > 1.0 / 86400.0) {
                // A well held at nothing lifts nothing: at q -> 0 the table
                // wants a full column, which a zero share never asks for.
                fail("tubing cannot lift the rate");     // rate/bhp/group control needs the tubing to allow more
                v.wells_unliftable.push_back(w);
            }
        }
        if (c == 'O' && well.oil_rate_limit > 0.0 && std::abs(qo - well.oil_rate_limit) > r_tol * well.oil_rate_limit) { fail("rate well off its limit"); }
        if (c == 'B' && std::abs(bhp[w] - well.bhp_limit) > dp_tol) { fail("bhp well off its limit"); }
        if (well.in_group) { in_group += well.efficiency * qo; any_held = any_held || c == 'G' || c == 'R'; }
    }
    if (sys.groupTarget() > 0.0) {
        if (in_group > sys.groupTarget() * (1 + r_tol)) { fail("group total above target"); }
        if (any_held && std::abs(in_group - sys.groupTarget()) > r_tol * sys.groupTarget()) { fail("held wells but target not met"); }
    }
    // The tree: no group above its target on its mode; a held well has a
    // bound ancestor -- one at its target -- and produces no more than
    // its own limits allow.
    if (sys.usesGroupTree()) {
        const auto& groups = sys.groups();
        std::vector<double> on_mode(sys.numGroups(), 0.0);
        for (int g = 0; g < sys.numGroups(); ++g) {
            const auto c = Sys::modeWeights(groups[g].mode, groups[g].resv_coeff);
            for (int w = 0; w < sys.numWells(); ++w) {
                // GEFAC applies where a group hands its rate to its parent,
                // as the group rows have it: q_g = sum_c eff_c q_c + sum_w eff_w q_w.
                bool under = false;
                double eff = sys.wells()[w].efficiency;
                for (int a = sys.wells()[w].group; a >= 0; a = groups[a].parent) {
                    if (a == g) { under = true; break; }
                    eff *= groups[a].efficiency;
                }
                if (!under) { continue; }
                for (int ph = 0; ph < NP; ++ph) { on_mode[g] += c[ph] * eff * q[w][ph]; }
            }
            // Relative, with a floor of 1 sm3/d: a zero target ("produce
            // nothing") is otherwise violated by rounding.
            if (groups[g].target > 0.0 && on_mode[g] > groups[g].target + std::max(r_tol * groups[g].target, 1.0 / 86400.0)) {
                // Over with every well under it on a control it cannot leave (not group controllable,
                // or no liftable share) is the limit's RATE action exhausted, not a wrong answer.
                bool reducible = false;
                for (int w = 0; w < sys.numWells() && !reducible; ++w) {
                    if (!(q_oil[w] > 0.0) || controls[w] == 'R' || controls[w] == 'G' || !sys.holdable(w)) { continue; }
                    for (int a = sys.wells()[w].group; a >= 0; a = groups[a].parent) {
                        if (a == g) { reducible = true; break; }
                    }
                }
                if (reducible) {
                    fail("group above its target");
                    v.groups_over.push_back(g);
                    v.over_ratio.push_back(on_mode[g] / groups[g].target);
                } else {
                    ++v.over_irreducible;
                }
            }
        }
        for (int w = 0; w < sys.numWells(); ++w) {
            if (controls[w] != 'R') { continue; }
            bool bound = false;
            for (int a = sys.wells()[w].group; a >= 0; a = groups[a].parent) {
                if (groups[a].target > 0.0 && std::abs(on_mode[a] - groups[a].target) <= std::max(r_tol * groups[a].target, 1.0 / 86400.0)) { bound = true; break; }
            }
            if (!bound) {
                std::string near;
                for (int a = sys.wells()[w].group; a >= 0; a = groups[a].parent) {
                    if (groups[a].target > 0.0) { near = fmt::format(" ({} at {:.2f} of target)", groups[a].name, on_mode[a] / groups[a].target); break; }
                }
                fail("held well without a group at its target" + near);
            }
            const auto& well = sys.wells()[w];
            if (well.vfp_table > 0) {
                const double cap = sys.thpPotential(well, p[well.node]);
                if (cap > 0.0 && cap < 1e30 && q_oil[w] > cap * (1 + r_tol)) { fail("held well above its tubing capacity"); }
            }
        }
    }
    // The network: pressures from the rates, top down.
    std::vector<std::array<double, NP>> Q(sys.numNodes() + 1);
    for (int n = sys.numNodes(); n >= 1; --n) {
        auto sum = sys.nodeSource(n);
        for (int w = 0; w < sys.numWells(); ++w) {
            if (sys.wells()[w].node != n) { continue; }
            for (int ph = 0; ph < NP; ++ph) { sum[ph] += sys.wells()[w].efficiency * (q[w][ph] + (ph == 2 ? sys.wells()[w].lift_gas : 0.0)); }
        }
        for (int c = 1; c <= sys.numNodes(); ++c) {
            if (sys.nodes()[c].parent != n) { continue; }
            for (int ph = 0; ph < NP; ++ph) { sum[ph] += sys.nodes()[c].efficiency * Q[c][ph]; }
        }
        Q[n] = sum;
    }
    for (int n = 1; n <= sys.numNodes(); ++n) {
        const auto& node = sys.nodes()[n];
        const double up = node.parent == 0 ? sys.terminalPressure() : p[node.parent];
        const double want = node.vfp_table != NetworkSolve::NoTable
            ? sys.tableBhp(node.vfp_table, up, Q[n], sys.branchAlq(n)) : up;
        if (std::abs(p[n] - want) > dp_tol) { fail("node pressure off its branch"); }
    }
    return v;
}

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_JUDGE_HEADER_INCLUDED
