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
#ifndef OPM_NETWORK_TUBING_EXTENSION_HEADER_INCLUDED
#define OPM_NETWORK_TUBING_EXTENSION_HEADER_INCLUDED

#include <opm/simulators/wells/network/NetworkReducedSolve.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace Opm::NetworkSolve {

template<class Scalar>
struct ExtensionResult
{
    bool converged = false;   // a solve converged with no open well on the continuation
    int passes = 0;
    int closed = 0;           // wells shut by the closing rule, all passes
    int iterations = 0;
    int evaluations = 0;
    ReducedResult<Scalar> last;
    std::vector<char> closed_wells;
    std::vector<Scalar> closed_gap;   // how far the line missed the curve, where each was closed
    std::vector<Scalar> closed_bhp;   // the bhp of the continuation point it was closed at
};

/// Which of the wells on the continuation a pass shuts. All: every one at
/// once. Sequential: the one whose line misses the curve by most, so the
/// others see the lower node pressure before they are judged. Tiered: every
/// well whose gap exceeds the table's own resolution at its touching vertex
/// at once, since the table is sure of those, then the marginal ones one at
/// a time. The decision is kept in one place so that reservoir or global
/// information can enter it later.
enum class Closing { All, Sequential, Tiered };

/// The reduced solve on the continued tubing curves (Stein's suggestion),
/// with the shut decision taken afterwards for every well at once: a well
/// whose converged point is on the continuation, where the real tubing has
/// no crossing, is shut, sticky, and the system solved again, until no open
/// well is on it. Shuts only, never revives, so it terminates. The
/// continuation is switched off again on return. keep_dead starts from the
/// system's committed dead set instead of clearing it.
///
/// Not a per-column change of the table data: the hump's flow moves with
/// the fraction columns (index 0 to 13 of 21 on the MODEL5 table), so a
/// flattened table interpolates to a different curve on the stable branch
/// too, and moved a flowing well's crossing by 1.3 %. Nor a running
/// minimum of the interpolated curve: its shallow dips above the hump are
/// stable against the line's slope, and flattening them moved a crossing
/// by 1 %. The continuation is anchored at the line's own touching point;
/// see ProductionSystem::setTubingExtension.
template<class Sys>
ExtensionResult<typename Sys::ScalarType>
solveReducedOnExtension(Sys& system,
                        const std::vector<typename Sys::ScalarType>& node_pressure_guess,
                        const Parameters<typename Sys::ScalarType> params,
                        const Closing closing = Closing::All,
                        const int max_passes = 20,
                        const bool keep_dead = false)
{
    using Scalar = typename Sys::ScalarType;
    constexpr int NP = Sys::NP;
    ExtensionResult<Scalar> out;
    out.closed_wells.assign(system.numWells(), 0);
    out.closed_gap.assign(system.numWells(), Scalar{0});
    out.closed_bhp.assign(system.numWells(), Scalar{0});
    system.setTubingExtension(true);
    auto p = node_pressure_guess;
    struct Candidate { int w; Scalar gap; Scalar res; Scalar bhp; };
    for (int pass = 0; pass < max_passes; ++pass) {
        auto r = solveReduced(system, p, params, true, CliffRule::Die, keep_dead || pass > 0);
        ++out.passes;
        out.iterations += r.iterations;
        out.evaluations += r.evaluations;
        out.last = std::move(r);
        if (!out.last.converged) { break; }
        p = out.last.node_pressure;
        std::vector<Candidate> on;
        const auto& dead = system.committedDead();
        for (int w = 0; w < system.numWells(); ++w) {
            const auto& well = system.wells()[w];
            if (well.shut || well.vfp_table <= 0 || system.control(w) == Sys::Control::Shut
                || (static_cast<std::size_t>(w) < dead.size() && dead[w])) { continue; }
            const Scalar qo = out.last.well_rate[w];
            if (!(qo > Scalar{0}) || !(well.ipr_b[1] < Scalar{0})) { continue; }
            const Scalar bhp = (qo - well.ipr_a[1]) / well.ipr_b[1];
            std::array<Scalar, NP> q{};
            for (int ph = 0; ph < NP; ++ph) { q[ph] = std::max(Sys::ipr(well, ph, bhp), Scalar{0}); }
            const Scalar p_node = well.node == 0 ? system.terminalPressure() : p[well.node];
            const auto tp = system.touchingPoint(well, p_node, q);
            if (tp.valid && tp.gap > Scalar{0}) { on.push_back({w, tp.gap, tp.res, bhp}); }
        }
        if (on.empty()) { out.converged = true; break; }
        std::sort(on.begin(), on.end(), [](const Candidate& a, const Candidate& b) { return a.gap > b.gap; });
        std::vector<Candidate> chosen;
        switch (closing) {
        case Closing::All:        chosen = on; break;
        case Closing::Sequential: chosen = {on.front()}; break;
        case Closing::Tiered:
            for (const auto& c : on) { if (c.gap > c.res) { chosen.push_back(c); } }
            if (chosen.empty()) { chosen = {on.front()}; }
            break;
        }
        for (const auto& c : chosen) {
            system.killWell(c.w);
            out.closed_wells[c.w] = 1;
            out.closed_gap[c.w] = c.gap;
            out.closed_bhp[c.w] = c.bhp;
            ++out.closed;
        }
    }
    system.setTubingExtension(false);
    return out;
}

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_TUBING_EXTENSION_HEADER_INCLUDED
