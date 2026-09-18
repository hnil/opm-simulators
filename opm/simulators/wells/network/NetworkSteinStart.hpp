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
#ifndef OPM_NETWORK_STEIN_START_HEADER_INCLUDED
#define OPM_NETWORK_STEIN_START_HEADER_INCLUDED

#include <opm/simulators/wells/network/NetworkProductionSystem.hpp>
#include <opm/simulators/wells/ProdGroupTreeBalancer.hpp>
#include <opm/simulators/wells/ProdGroupTreeNode.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>

#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <algorithm>
#include <array>
#include <string>
#include <vector>

namespace Opm::NetworkSolve {

/// Stein's balancer on the system's own tree, with each well's capacity taken
/// at the state `x` -- the same IPR, tubing and pressures the Newton sees.
/// `guide_rate` is filled here from the system's guides and must not be the
/// simulator's own object. Oil rates per well (SI); empty if the balancer
/// rejected the tree. Wells outside the tree are left at their allowance.
template<class Scalar>
std::vector<Scalar> steinAllocation(ProductionSystem<Scalar>& sys,
                                    const std::vector<Scalar>& x,
                                    GuideRate& guide_rate,
                                    const int step)
{
    using Sys = ProductionSystem<Scalar>;
    sys.updateControls(x);                 // the capacities at these pressures
    DeferredLogger logger;
    const auto& groups = sys.groups();
    if (groups.empty()) {
        return {};
    }
    // The guide as the oil potential and the other phases in the well's own
    // proportions at its bhp limit, so a group's guide on any mode is the sum
    // of its wells' -- the system's rule, read back through GuideRate.
    std::vector<std::array<double, 3>> gpot(sys.numGroups(), {0.0, 0.0, 0.0});   // oil, gas, water
    for (int w = 0; w < sys.numWells(); ++w) {
        const auto& well = sys.wells()[w];
        if (well.group < 0) {
            continue;
        }
        const double bhp = well.bhp_limit;
        const double oil = std::max(double(well.ipr_a[1] + well.ipr_b[1] * bhp), 1e-30);
        const double guide = std::max(double(well.guide), 1e-9 / 86400.0);
        const std::array<double, 3> pot{guide,
                                        guide * std::max(double(well.ipr_a[2] + well.ipr_b[2] * bhp), 0.0) / oil,
                                        guide * std::max(double(well.ipr_a[0] + well.ipr_b[0] * bhp), 0.0) / oil};
        guide_rate.compute(well.name, step, 0.0, pot[0], pot[1], pot[2]);
        for (int g = well.group; g >= 0; g = groups[g].parent) {
            for (int i = 0; i < 3; ++i) { gpot[g][i] += well.efficiency * pot[i]; }
        }
    }
    for (int g = 0; g < sys.numGroups(); ++g) {
        guide_rate.compute(groups[g].name, step, 0.0, gpot[g][0], gpot[g][1], gpot[g][2]);
    }
    auto wmode = [](const typename Sys::Mode m) {
        switch (m) {
        case Sys::Mode::Water:  return ::Opm::Well::ProducerCMode::WRAT;
        case Sys::Mode::Gas:    return ::Opm::Well::ProducerCMode::GRAT;
        case Sys::Mode::Liquid: return ::Opm::Well::ProducerCMode::LRAT;
        default:                return ::Opm::Well::ProducerCMode::ORAT;
        }
    };
    auto gmode = [](const typename Sys::Mode m) {
        switch (m) {
        case Sys::Mode::Water:  return ::Opm::Group::ProductionCMode::WRAT;
        case Sys::Mode::Gas:    return ::Opm::Group::ProductionCMode::GRAT;
        case Sys::Mode::Liquid: return ::Opm::Group::ProductionCMode::LRAT;
        default:                return ::Opm::Group::ProductionCMode::ORAT;
        }
    };
    ProdGroupTreeBalancer::Tree<Scalar> tree;
    for (int g = 0; g < sys.numGroups(); ++g) {
        ProdGroupTreeNode<Scalar> n;
        n.name = groups[g].name;
        n.type = ProdNodeType::Group;
        n.parent = groups[g].parent >= 0 ? groups[groups[g].parent].name : std::string{};
        n.availableForGroupControl = groups[g].available;
        n.efficiencyFactor = groups[g].efficiency;
        n.hasGuideRate = false;
        const bool has_target = groups[g].target > Scalar{0};
        n.modeCategory = has_target ? ProdNodeModeCategory::Individual : ProdNodeModeCategory::Group;
        if (has_target) {
            n.mode = wmode(groups[g].mode);
            n.preferredMode = gmode(groups[g].mode);
            n.Limits[n.mode] = groups[g].target;
        }
        tree.emplace(n.name, std::move(n));
    }
    for (int g = 0; g < sys.numGroups(); ++g) {
        if (groups[g].parent >= 0) {
            tree.at(groups[groups[g].parent].name).children.push_back(groups[g].name);
        }
    }
    for (int w = 0; w < sys.numWells(); ++w) {
        const auto& well = sys.wells()[w];
        if (well.group < 0) {
            continue;
        }
        ProdGroupTreeNode<Scalar> n;
        n.name = well.name;
        n.type = ProdNodeType::Well;
        n.parent = groups[well.group].name;
        n.availableForGroupControl = true;
        n.efficiencyFactor = well.efficiency;
        n.mode = ::Opm::Well::ProducerCMode::GRUP;
        n.hasGuideRate = true;
        const Scalar allow = well.shut ? Scalar{0} : std::max(sys.ownAllowance(w), Scalar{0});
        n.Limits[::Opm::Well::ProducerCMode::ORAT] = allow;
        std::array<Scalar, Sys::NP> q{};           // water, oil, gas
        if (well.ipr_b[1] < Scalar{0}) {
            const Scalar bhp = (allow - well.ipr_a[1]) / well.ipr_b[1];
            for (int ph = 0; ph < Sys::NP; ++ph) {
                q[ph] = std::max(well.ipr_a[ph] + well.ipr_b[ph] * bhp, Scalar{0});
            }
        } else {
            q[1] = allow;
        }
        n.rates = {-q[1], -q[0], -q[2]};           // canonical [oil, water, gas], production negative
        n.initialRates = n.rates;
        tree.at(n.parent).children.push_back(n.name);
        tree.emplace(well.name, std::move(n));
    }
    if (!ProdGroupTreeBalancer::balanceTreeForTesting(tree, guide_rate, Scalar{1e-8}, logger)) {
        return {};
    }
    std::vector<Scalar> q_oil(sys.numWells());
    for (int w = 0; w < sys.numWells(); ++w) {
        const auto& well = sys.wells()[w];
        q_oil[w] = (well.group < 0) ? std::max(sys.ownAllowance(w), Scalar{0})
                                    : -tree.at(well.name).rates[0];
    }
    return q_oil;
}

template<class Scalar>
struct SteinStart
{
    bool ok = false;
    std::vector<Scalar> node_pressure;     // index 0 the terminal
    std::vector<Scalar> well_rate;         // oil, SI
    int sweeps = 0;
};

/// The starting point for a solve with the IPR as given: capacities at the
/// guessed pressures, Stein's allocation from them, and the node pressures
/// those rates imply. Repeated `sweeps` times -- each sweep is one step of the
/// simulator's own fixed point, with Stein's balancer as the well model. The
/// caller decides whether the rates also open the wells (setStartAllocation):
/// the full routes start from them, the reduced route only from the pressures.
template<class Scalar>
SteinStart<Scalar> steinStart(ProductionSystem<Scalar>& sys,
                              const std::vector<Scalar>& node_pressure_guess,
                              GuideRate& guide_rate,
                              const int step,
                              const int sweeps = 1)
{
    SteinStart<Scalar> out;
    out.node_pressure = node_pressure_guess;
    for (int k = 0; k < sweeps; ++k) {
        const auto x = sys.start(out.node_pressure);
        auto q = steinAllocation(sys, x, guide_rate, step);
        if (q.empty()) {
            return out;
        }
        out.node_pressure = sys.pressuresFromWellRates(q);
        out.well_rate = std::move(q);
        out.ok = true;
        out.sweeps = k + 1;
    }
    return out;
}

} // namespace Opm::NetworkSolve

#endif // OPM_NETWORK_STEIN_START_HEADER_INCLUDED
