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

#include <config.h>

#define BOOST_TEST_MODULE ProdGroupTreeBalancerTest

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/ProdGroupTreeBalancer.hpp>
#include <opm/simulators/wells/ProdGroupTreeNode.hpp>

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

// Trees built by hand, in sm3/d, balanced through balanceTreeForTesting. The
// schedule exists only to give GuideRate somewhere to keep the potentials.

namespace {

using Tree = Opm::ProdGroupTreeBalancer::Tree<double>;
using Node = Opm::ProdGroupTreeNode<double>;
using Opm::Well;
using Opm::Group;

struct Fixture
{
    std::unique_ptr<Opm::Deck> deck;
    std::unique_ptr<Opm::EclipseState> es;
    std::unique_ptr<Opm::Schedule> schedule;
    std::unique_ptr<Opm::GuideRate> guide;
    Tree tree;

    explicit Fixture(const std::string& sched)
    {
        std::string d;
        d += "RUNSPEC\nTITLE\nBALANCER_UNIT\nDIMENS\n 20 1 2 /\nOIL\nWATER\nGAS\nDISGAS\nMETRIC\n";
        d += "START\n 1 'JAN' 2020 /\nWELLDIMS\n 4 2 6 4 /\nEQLDIMS\n 1 1* 25 1* 1 /\n";
        d += "TABDIMS\n 1 1 50 60 1 60 1 1 /\nUNIFOUT\n";
        d += "GRID\nINIT\nDXV\n20*1000 /\nDYV\n1*1000 /\nDZV\n50 50 /\nTOPS\n20*7000 /\n";
        d += "PORO\n40*0.2 /\nPERMX\n40*200 /\nPERMY\n40*200 /\nPERMZ\n40*20 /\n";
        d += "PROPS\nSWOF\n 0.2 0 1 0\n 0.3 0.07 0.8 0\n 1.0 1 0 0 /\nSGOF\n 0 0 1 0\n 0.05 0 0.8 0\n 0.79 1 0 0 /\n";
        d += "DENSITY\n 800 1000 1 /\nPVTW\n 1 1.0 4.0E-5 0.5 0.0 /\nPVDG\n 1 1.0 0.01\n 100 0.1 0.015\n 300 0.033 0.02 /\n";
        d += "PVTO\n 1 50 1.2 1.0\n 150 1.15 1.1\n 300 1.10 1.2 /\n 10 150 1.25 0.9\n 250 1.20 1.0\n 350 1.15 1.1 /\n/\n";
        d += "REGIONS\nEQLNUM\n40*1 /\nSOLUTION\nEQUIL\n 7000 270 7050 0 7000 0 1* 0 0 /\nSUMMARY\nFOPR\nSCHEDULE\n";
        d += sched + "TSTEP\n1 /\nEND\n";
        deck = std::make_unique<Opm::Deck>(Opm::Parser{}.parseString(d));
        es = std::make_unique<Opm::EclipseState>(*deck);
        schedule = std::make_unique<Opm::Schedule>(*deck, *es);
        guide = std::make_unique<Opm::GuideRate>(*schedule);
        Node field;
        field.name = "FIELD"; field.type = Opm::ProdNodeType::Group;
        field.availableForGroupControl = false; field.modeCategory = Opm::ProdNodeModeCategory::Group;
        tree.emplace("FIELD", field);
    }

    void group(const std::string& name, const std::string& parent,
               Well::ProducerCMode mode = Well::ProducerCMode::CMODE_UNDEFINED, double target = 0.0)
    {
        Node n;
        n.name = name; n.type = Opm::ProdNodeType::Group; n.parent = parent;
        n.availableForGroupControl = true; n.hasGuideRate = false;
        if (target > 0.0) {
            n.modeCategory = Opm::ProdNodeModeCategory::Individual;
            n.mode = mode;
            n.preferredMode = mode == Well::ProducerCMode::GRAT ? Group::ProductionCMode::GRAT
                                                                : Group::ProductionCMode::ORAT;
            n.Limits[mode] = target;
        } else {
            n.modeCategory = Opm::ProdNodeModeCategory::Group;
        }
        tree.at(parent).children.push_back(name);
        tree.emplace(name, n);
    }

    // A well with its current rates {oil, water, gas} and its own limits.
    void well(const std::string& name, const std::string& parent, std::array<double, 3> rates,
              std::map<Well::ProducerCMode, double> limits)
    {
        Node n;
        n.name = name; n.type = Opm::ProdNodeType::Well; n.parent = parent;
        n.availableForGroupControl = true; n.hasGuideRate = true;
        n.mode = Well::ProducerCMode::GRUP;
        n.Limits = std::move(limits);
        n.rates = {-rates[0], -rates[1], -rates[2]};
        n.initialRates = n.rates;
        tree.at(parent).children.push_back(name);
        tree.emplace(name, n);
        // Equal guides: the potential is the current rate.
        guide->compute(name, 0, 0.0, rates[0], rates[2], rates[1]);
    }

    bool balance(const bool assignTargets = false)
    {
        // Groups get the sum of their wells as potential, bottom-up.
        for (const auto& [name, node] : tree) {
            if (node.type != Opm::ProdNodeType::Group) continue;
            std::array<double, 3> pot{0.0, 0.0, 0.0};
            std::vector<std::string> stack = node.children;
            while (!stack.empty()) {
                const auto c = stack.back(); stack.pop_back();
                const auto& cn = tree.at(c);
                if (cn.type == Opm::ProdNodeType::Well) {
                    for (int i = 0; i < 3; ++i) pot[i] += -cn.rates[i];
                } else {
                    stack.insert(stack.end(), cn.children.begin(), cn.children.end());
                }
            }
            guide->compute(name, 0, 0.0, pot[0], pot[2], pot[1]);
        }
        Opm::DeferredLogger logger;
        return Opm::ProdGroupTreeBalancer::balanceTreeForTesting(tree, *guide, 1e-8, logger, assignTargets);
    }

    double oil(const std::string& name) const { return -tree.at(name).rates[0]; }
    double gas(const std::string& name) const { return -tree.at(name).rates[2]; }
};

const std::string kThreeWells =
    "GRUPTREE\n 'PLAT' 'FIELD' /\n 'G1' 'PLAT' /\n 'G3' 'PLAT' /\n 'G5' 'G3' /\n/\n"
    "WELSPECS\n 'W1' 'G1' 1 1 7000 'OIL' /\n 'W2' 'G5' 2 1 7000 'OIL' /\n 'W3' 'G5' 3 1 7000 'OIL' /\n/\n"
    "COMPDAT\n 'W1' 1 1 1 2 'OPEN' 1* 1* 0.2 /\n 'W2' 2 1 1 2 'OPEN' 1* 1* 0.2 /\n 'W3' 3 1 1 2 'OPEN' 1* 1* 0.2 /\n/\n"
    "WCONPROD\n 'W1' 'OPEN' 'GRUP' 2000 4* 100 /\n 'W2' 'OPEN' 'GRUP' 2000 4* 100 /\n 'W3' 'OPEN' 'GRUP' 2000 4* 100 /\n/\n"
    "GCONPROD\n 'PLAT' 'ORAT' 3000 3* 'RATE' /\n 'G5' 'ORAT' 500 3* 'RATE' /\n/\n";

} // anonymous namespace

// A group's own target two levels under a bound parent, with a group without a
// guide rate between: G5 may make 500 of PLAT's 3000, W1 2000, so PLAT is short
// at 2500. The balancer used to hand G5 2000, because with no current rate on
// G5 its own limit was never applied to the target it was given.
BOOST_AUTO_TEST_CASE(a_nested_target_under_a_transparent_group)
{
    Fixture f(kThreeWells);
    f.group("PLAT", "FIELD", Well::ProducerCMode::ORAT, 3000.0);
    f.group("G1", "PLAT");
    f.group("G3", "PLAT");
    f.group("G5", "G3", Well::ProducerCMode::ORAT, 500.0);
    for (const auto& [w, g] : {std::pair{"W1", "G1"}, std::pair{"W2", "G5"}, std::pair{"W3", "G5"}}) {
        f.well(w, g, {2000.0, 0.0, 200000.0}, {{Well::ProducerCMode::ORAT, 2000.0}});
    }
    BOOST_REQUIRE(f.balance());
    BOOST_CHECK_CLOSE(f.oil("W2") + f.oil("W3"), 500.0, 1e-6);
    BOOST_CHECK_CLOSE(f.oil("W1"), 2000.0, 1e-6);
    BOOST_CHECK_CLOSE(f.oil("G5"), 500.0, 1e-6);
    BOOST_CHECK(f.tree.at("G5").modeCategory == Opm::ProdNodeModeCategory::Individual);
}

// A well carrying two of its own limits: its share of the group target is
// feasible on oil but not on gas, so it is held at its GRAT limit and the
// sibling takes the rest. This is what the tree has to see for a well's
// second limit to be decided at allocation time rather than by the well solve.
BOOST_AUTO_TEST_CASE(a_well_with_two_limits_is_held_at_the_binding_one)
{
    Fixture f("GRUPTREE\n 'PLAT' 'FIELD' /\n 'G1' 'PLAT' /\n/\n"
              "WELSPECS\n 'W1' 'G1' 1 1 7000 'OIL' /\n 'W2' 'G1' 2 1 7000 'OIL' /\n/\n"
              "COMPDAT\n 'W1' 1 1 1 2 'OPEN' 1* 1* 0.2 /\n 'W2' 2 1 1 2 'OPEN' 1* 1* 0.2 /\n/\n"
              "WCONPROD\n 'W1' 'OPEN' 'GRUP' 2000 4* 100 /\n 'W2' 'OPEN' 'GRUP' 2000 1* 250000 2* 100 /\n/\n"
              "GCONPROD\n 'PLAT' 'ORAT' 3000 3* 'RATE' /\n/\n");
    f.group("PLAT", "FIELD", Well::ProducerCMode::ORAT, 3000.0);
    f.group("G1", "PLAT");
    f.well("W1", "G1", {1500.0, 0.0, 150000.0}, {{Well::ProducerCMode::ORAT, 2000.0}});
    f.well("W2", "G1", {1500.0, 0.0, 300000.0},
           {{Well::ProducerCMode::ORAT, 2000.0}, {Well::ProducerCMode::GRAT, 250000.0}});
    BOOST_REQUIRE(f.balance());
    BOOST_CHECK_CLOSE(f.gas("W2"), 250000.0, 1e-6);
    BOOST_CHECK_CLOSE(f.oil("W2"), 1250.0, 1e-6);
    BOOST_CHECK_CLOSE(f.oil("W1"), 1750.0, 1e-6);
    BOOST_CHECK_CLOSE(f.oil("PLAT"), 3000.0, 1e-6);
    BOOST_CHECK(f.tree.at("W2").modeCategory == Opm::ProdNodeModeCategory::Individual);
    BOOST_CHECK(f.tree.at("W2").mode == Well::ProducerCMode::GRAT);
    BOOST_CHECK(f.tree.at("W1").modeCategory == Opm::ProdNodeModeCategory::Group);
}

// A well that cannot be held already makes more than its group's target: nothing is
// left to share, and the held well gets zero -- not a negative rate (GRPFLD-02:
// -6.2e6 sm3/d on B-2H, and the route read the well as not held at all).
BOOST_AUTO_TEST_CASE(nothing_left_to_share_gives_zero_not_a_negative_rate)
{
    Fixture f("GRUPTREE\n 'PLAT' 'FIELD' /\n 'G1' 'PLAT' /\n/\n"
              "WELSPECS\n 'W1' 'G1' 1 1 7000 'OIL' /\n 'W2' 'G1' 2 1 7000 'OIL' /\n/\n"
              "COMPDAT\n 'W1' 1 1 1 2 'OPEN' 1* 1* 0.2 /\n 'W2' 2 1 1 2 'OPEN' 1* 1* 0.2 /\n/\n"
              "WCONPROD\n 'W1' 'OPEN' 'ORAT' 3500 4* 100 /\n 'W2' 'OPEN' 'GRUP' 2000 4* 100 /\n/\n"
              "GCONPROD\n 'PLAT' 'ORAT' 3000 3* 'RATE' /\n/\n");
    f.group("PLAT", "FIELD", Well::ProducerCMode::ORAT, 3000.0);
    f.group("G1", "PLAT");
    f.well("W1", "G1", {3500.0, 0.0, 350000.0}, {{Well::ProducerCMode::ORAT, 3500.0}});
    auto& w1 = f.tree.at("W1");
    w1.availableForGroupControl = false;
    w1.modeCategory = Opm::ProdNodeModeCategory::Individual;
    w1.mode = Well::ProducerCMode::ORAT;
    f.well("W2", "G1", {1000.0, 0.0, 100000.0}, {{Well::ProducerCMode::ORAT, 2000.0}});
    f.balance();
    BOOST_TEST_MESSAGE("W2 oil " << f.oil("W2") << " category " << static_cast<int>(f.tree.at("W2").modeCategory));
    BOOST_CHECK_GE(f.oil("W2"), 0.0);
    BOOST_CHECK_SMALL(f.oil("W2"), 1e-6);
    BOOST_CHECK_CLOSE(f.oil("W1"), 3500.0, 1e-6);
    // Still holding its limit, with nothing left to give: the held well hangs off it.
    BOOST_CHECK(f.tree.at("PLAT").modeCategory == Opm::ProdNodeModeCategory::Individual);
    BOOST_CHECK(f.tree.at("W2").modeCategory == Opm::ProdNodeModeCategory::Group);
}

// With assignTargets the allocation is written as each well's group target: a
// group-controlled well's target is its allocated rate in the group's mode.
BOOST_AUTO_TEST_CASE(assign_targets_writes_the_allocation)
{
    Fixture f(kThreeWells);
    f.group("PLAT", "FIELD", Well::ProducerCMode::ORAT, 3000.0);
    f.group("G1", "PLAT");
    f.group("G3", "PLAT");
    f.group("G5", "G3", Well::ProducerCMode::ORAT, 500.0);
    for (const auto& [w, g] : {std::pair{"W1", "G1"}, std::pair{"W2", "G5"}, std::pair{"W3", "G5"}}) {
        f.well(w, g, {2000.0, 0.0, 200000.0}, {{Well::ProducerCMode::ORAT, 2000.0}});
    }
    BOOST_REQUIRE(f.balance(/*assignTargets*/ true));
    for (const auto& w : {"W2", "W3"}) {
        const auto& gt = f.tree.at(w).groupTarget;
        BOOST_CHECK(f.tree.at(w).modeCategory == Opm::ProdNodeModeCategory::Group);
        BOOST_CHECK(gt.ctrlMode == Group::ProductionCMode::ORAT);
        BOOST_CHECK_EQUAL(gt.groupName, "G5");
        BOOST_CHECK_CLOSE(gt.value, 250.0, 1e-6);
    }
    // W1 is at its own limit, and PLAT ends short of its target with no
    // group-controlled child, so W1 is under no group control at all.
    BOOST_CHECK(f.tree.at("W1").modeCategory == Opm::ProdNodeModeCategory::Individual);
    BOOST_CHECK(f.tree.at("PLAT").modeCategory == Opm::ProdNodeModeCategory::None);
    BOOST_CHECK(f.tree.at("W1").groupTarget.ctrlMode == Group::ProductionCMode::NONE);
}

// A well with no gas under a gas target: the target cannot hold it, so it stays
// where its own limits put it and the gas producer beside it carries the target.
// The balancer used to take the gas target as the dry well's own and zero it.
BOOST_AUTO_TEST_CASE(a_well_without_the_target_phase_is_not_zeroed)
{
    Fixture f("GRUPTREE\n 'PLAT' 'FIELD' /\n 'G1' 'PLAT' /\n 'G2' 'PLAT' /\n 'M2' 'G2' /\n/\n"
              "WELSPECS\n 'W1' 'G1' 1 1 7000 'OIL' /\n 'W2' 'M2' 2 1 7000 'OIL' /\n/\n"
              "COMPDAT\n 'W1' 1 1 1 2 'OPEN' 1* 1* 0.2 /\n 'W2' 2 1 1 2 'OPEN' 1* 1* 0.2 /\n/\n"
              "WCONPROD\n 'W1' 'OPEN' 'GRUP' 2000 4* 100 /\n 'W2' 'OPEN' 'GRUP' 1000 4* 100 /\n/\n"
              "GCONPROD\n 'PLAT' 'GRAT' 2* 200000 1* 'RATE' /\n/\n");
    f.group("PLAT", "FIELD", Well::ProducerCMode::GRAT, 200000.0);
    f.group("G1", "PLAT");
    // G2 answers to its own limits only (GCONPROD item 8 NO), gas the preferred one.
    f.group("G2", "PLAT", Well::ProducerCMode::GRAT, 800000.0);
    f.tree.at("G2").availableForGroupControl = false;
    f.tree.at("G2").Limits[Well::ProducerCMode::ORAT] = 5000.0;
    // As the tree builder leaves a group: no category of its own yet.
    f.tree.at("G2").modeCategory = Opm::ProdNodeModeCategory::Group;
    f.tree.at("G2").mode = Well::ProducerCMode::CMODE_UNDEFINED;
    f.well("W1", "G1", {1500.0, 0.0, 300000.0}, {{Well::ProducerCMode::ORAT, 2000.0}});
    // The dry well as the simulator hands it over: on its own control, at a capacity
    // given as a limit on the sum of its rates, with a gas limit of its own as well.
    f.group("M2", "G2");     // no limits, no guide rate: between the group and its well
    // Its guide rate comes from potentials that do have gas; the rates it is handed
    // over with do not.
    f.well("W2", "M2", {1000.0, 0.5, 200000.0},
           {{Well::ProducerCMode::ORAT, 4000.0}, {Well::ProducerCMode::GRAT, 500000.0},
            {Well::ProducerCMode::LRAT, 8000.0}, {Well::ProducerCMode::THP, 1000.5}});
    f.tree.at("W2").modeCategory = Opm::ProdNodeModeCategory::Individual;
    f.tree.at("W2").mode = Well::ProducerCMode::THP;
    f.tree.at("W2").rates = {-1000.0, -0.5, 0.0};
    f.tree.at("W2").initialRates = f.tree.at("W2").rates;
    BOOST_REQUIRE(f.balance());
    BOOST_CHECK_CLOSE(f.oil("W2"), 1000.0, 1e-6);
    BOOST_CHECK_CLOSE(f.gas("W1"), 200000.0, 1e-6);
    BOOST_CHECK_CLOSE(f.oil("W1"), 1000.0, 1e-6);
}
