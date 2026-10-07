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

#include <opm/simulators/wells/network/NetworkSolve.hpp>

#include <cstddef>
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

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_GAS_LIFT_HEADER_INCLUDED
