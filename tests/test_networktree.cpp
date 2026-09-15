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
#include "NetworkSolveTestSupport.hpp"

// The group tree and the reduced form: the tree as rows, against Stein's
// balancer, and route C against the outer loop it is the Schur complement of.

BOOST_AUTO_TEST_SUITE(NetworkTreeBench)

// The sparse group equations against Stein's ProdGroupTreeBalancer, on the
// same deck.
//
// The two reach the allocation by different routes -- his by a sorted
// guide-rate-to-limit recursion, this by complementarity with multipliers --
// under the same fixed-phase-fraction assumption. So where his accepts a
// configuration, the equations should produce the same rates, and a
// disagreement is informative rather than ambiguous.
BOOST_AUTO_TEST_CASE(the_group_tree_against_steins_balancer)
{
    const auto deck_path = std::filesystem::path(__FILE__).parent_path() / "GROUPTREE.DATA";
    if (!std::filesystem::exists(deck_path)) {
        BOOST_TEST_MESSAGE("GROUPTREE.DATA not found, skipping");
        return;
    }
    Parser parser;
    const auto deck = parser.parseFile(deck_path.string());
    const EclipseState es{deck};
    const Schedule schedule{deck, es};
    SummaryState summary_state{TimeService::now(), 0.0};

    const std::map<std::string, std::string> parent{
        {"PLAT", "FIELD"}, {"GP1", "PLAT"}, {"GP2", "PLAT"},
        {"W1", "GP1"}, {"W2", "GP1"}, {"W3", "GP2"}};
    const std::vector<std::string> wnames{"W1", "W2", "W3"};

    struct Case {
        const char* what;
        double target;                      // PLAT ORAT, sm3/d
        std::array<double, 3> cap;          // each well's own ORAT limit
        std::array<double, 3> guide;        // well guide rates
    };
    const std::vector<Case> cases{
        {"target binds, one well capped",   3000.0, {2000, 2000,  400}, {1, 1, 1}},
        {"target binds, nothing else",      1200.0, {2000, 2000, 2000}, {1, 1, 1}},
        {"target above what the wells make", 9000.0, {2000, 2000,  400}, {1, 1, 1}},
        {"two wells capped",                3000.0, {2000,  900,  400}, {1, 1, 1}},
        {"unequal guide rates",             3000.0, {2000, 2000, 2000}, {3, 2, 1}},
    };

    // Both routes to the active set are held to the same oracle: complementarity
    // inside the Newton, and the tree pass outside it.
    for (const bool active_set : {false, true}) {
    for (const auto& c : cases) {
        // ---- the oracle ---------------------------------------------------
        GuideRate guide_rate{schedule};
        DeferredLogger logger;
        // A fresh GuideRate holds nothing; the simulator fills it from the well
        // potentials each step. Without it every guide rate reads zero, the
        // sorted recursion has nothing to distribute by, and the balancer
        // allocates nothing while still reporting the tree valid.
        for (std::size_t i = 0; i < wnames.size(); ++i) {
            guide_rate.compute(wnames[i], /*report_step=*/0, /*sim_time=*/0.0,
                               c.guide[i] / 86400.0, 0.0, 0.0);
        }

        ProdGroupTreeBalancer::Tree<double> tree;
        auto addGroup = [&](const std::string& name, const std::vector<std::string>& kids,
                            const bool has_target) {
            ProdGroupTreeNode<double> n;
            n.name = name;
            n.type = ProdNodeType::Group;
            n.parent = parent.count(name) ? parent.at(name) : std::string{};
            n.children = kids;
            n.availableForGroupControl = true;
            n.modeCategory = has_target ? ProdNodeModeCategory::Individual
                                        : ProdNodeModeCategory::Group;
            if (has_target) {
                n.mode = Opm::Well::ProducerCMode::ORAT;
                n.preferredMode = Opm::Group::ProductionCMode::ORAT;
                n.Limits[Opm::Well::ProducerCMode::ORAT] = c.target / 86400.0;
            }
            tree.emplace(name, std::move(n));
        };
        addGroup("FIELD", {"PLAT"}, false);
        addGroup("PLAT", {"GP1", "GP2"}, true);
        addGroup("GP1", {"W1", "W2"}, false);
        addGroup("GP2", {"W3"}, false);
        for (std::size_t i = 0; i < wnames.size(); ++i) {
            ProdGroupTreeNode<double> n;
            n.name = wnames[i];
            n.type = ProdNodeType::Well;
            n.parent = parent.at(wnames[i]);
            n.availableForGroupControl = true;
            n.mode = Opm::Well::ProducerCMode::GRUP;
            // populateWellNode() sets this for every well; the propagation pass
            // only reads it, so without it no guide rate is ever assigned and
            // the balancer allocates nothing while calling the tree valid.
            n.hasGuideRate = true;
            n.Limits[Opm::Well::ProducerCMode::ORAT] = c.cap[i] / 86400.0;
            n.rates = {-c.cap[i] / 86400.0, 0.0, 0.0};   // [oil, water, gas], production negative
            n.initialRates = n.rates;
            tree.emplace(wnames[i], std::move(n));
        }
        const bool oracle_ok =
            ProdGroupTreeBalancer::balanceTreeForTesting(tree, guide_rate, 1e-8, logger);

        // ---- the equations --------------------------------------------------
        using Sys = NetworkSolve::ProductionSystem<double>;
        VFPProdProperties<double> props;
        const UnitSystem units{};
        Sys sys(props, units);
        NetworkSolve::Node term; term.name = "TERM"; term.parent = -1; sys.addNode(term, 0.0);
        NetworkSolve::Node nd; nd.name = "N"; nd.parent = 0; sys.addNode(nd, 0.0);
        sys.setTerminalPressure(convert::from(120.0, bars));

        // A group's guide rate is the sum of its children's. Leaving them at
        // the default 1 makes PLAT split evenly between GP1 (two wells) and
        // GP2 (one), which is a different answer from splitting evenly between
        // the three wells -- 300/300/600 against 400/400/400. The oracle
        // derives group guides from what is beneath; this has to match.
        const std::map<std::string, double> gguide{
            {"GP1", c.guide[0] + c.guide[1]}, {"GP2", c.guide[2]},
            {"PLAT", c.guide[0] + c.guide[1] + c.guide[2]},
            {"FIELD", c.guide[0] + c.guide[1] + c.guide[2]}};
        std::map<std::string, int> gidx;
        for (const auto& name : {"FIELD", "PLAT", "GP1", "GP2"}) {
            typename Sys::Group g;
            g.name = name;
            g.parent = parent.count(name) ? gidx.at(parent.at(name)) : -1;
            g.guide = gguide.at(name);
            if (std::string(name) == "PLAT") {
                g.mode = Sys::Mode::Oil;
                g.target = c.target / 86400.0;
            }
            gidx[name] = sys.addGroup(std::move(g));
        }
        for (std::size_t i = 0; i < wnames.size(); ++i) {
            typename Sys::Well w;
            w.name = wnames[i];
            w.node = 1;
            w.vfp_table = 0;
            const double cap = c.cap[i] / 86400.0;
            const double b = -cap / (250e5 - 100e5);
            w.ipr_a[1] = -b * 250e5;
            w.ipr_b[1] = b;
            w.bhp_limit = convert::from(100.0, bars);
            w.oil_rate_limit = cap;
            w.guide = c.guide[i];
            w.group = gidx.at(parent.at(wnames[i]));
            w.q_start = 0.5 * cap;
            sys.addWell(std::move(w));
        }
        sys.setGroupTree(true);
        sys.setAnalyticJacobian(true);
        sys.setComplementarity(true);
        sys.setGroupActiveSet(active_set);
        sys.finish();
        sys.finishGroups();

        const auto r = NetworkSolve::solve(sys, {convert::from(120.0, bars)},
                                           {1e-2, 80}, NetworkSolve::FullStep{});

        // ---- compare ---------------------------------------------------------
        std::string both;
        double worst = 0.0;
        for (int w = 0; w < sys.numWells(); ++w) {
            const auto& name = sys.wells()[w].name;
            const double mine = r.well_rate[w] * 86400.0;
            const double theirs = -tree.at(name).rates[0] * 86400.0;
            both += fmt::format(" {}={:.1f}/{:.1f}", name, mine, theirs);
            worst = std::max(worst, std::abs(mine - theirs) / std::max(std::abs(theirs), 1.0));
        }
        BOOST_TEST_MESSAGE((active_set ? "[active set] " : "[fb] ")
                           << c.what << ": equations/oracle" << both
                           << "  (" << r.iterations << " it, worst gap "
                           << 100 * worst << " %)");
        BOOST_CHECK(oracle_ok);
        BOOST_CHECK(r.converged);
        BOOST_CHECK_LT(worst, 0.01);
    }
    }
}

// The group rows' analytic derivatives against differences, entry by entry.
//
// This is the check that says the linearisation is right, independently of
// whether any solve converges: if the assembled Jacobian matches the
// difference quotient everywhere, the rows are differentiated correctly.
BOOST_AUTO_TEST_CASE(the_group_jacobian_matches_differences)
{
    TreeCase tc;
    struct Shape { const char* what; double field; double p1;
                   std::array<double,3> caps; TreeCase::Sys::Mode mode; double wct; };
    const std::vector<Shape> shapes{
        {"oil target",            3000.0,    0.0, {2000, 2000, 400}, TreeCase::Sys::Mode::Oil,    0.0},
        {"target at both levels", 3000.0, 2000.0, {2000, 2000, 400}, TreeCase::Sys::Mode::Oil,    0.0},
        {"liquid target",         3750.0,    0.0, {2000, 2000, 400}, TreeCase::Sys::Mode::Liquid, 0.25},
    };
    for (const auto& sh : shapes) {
        auto sys = tc.build(sh.field, sh.p1, sh.caps, sh.mode, sh.wct);
        auto x = sys.start({convert::from(120.0, bars)});
        // Away from the start, where the slacks are not all at their initial
        // values and the kinks are not sitting exactly on zero.
        for (int w = 0; w < sys.numWells(); ++w) {
            x[sys.wellRateIndex(w, 1)] *= 0.83;
            x[sys.wellBhpIndex(w)] += convert::from(4.0, bars);
        }
        // The group unknowns too: start() seeds the group rates from the
        // wells and parks every multiplier at the ceiling, which is a kink of
        // the ceiling's phi.
        for (int g = 0; g < sys.numGroups(); ++g) {
            for (int ph = 0; ph < TreeCase::Sys::NP; ++ph) { x[sys.gqIdx(g, ph)] *= 0.91; }
            x[sys.glIdx(g)] *= 0.7;
        }
        const auto J = sys.jacobian(x);
        const auto r0 = sys.residual(x);
        const int n = sys.size();
        double worst = 0.0;
        int worst_i = -1, worst_j = -1;
        for (int j = 0; j < n; ++j) {
            auto shifted = x;
            const double h = 1e-6 * sys.columnScale(j);
            shifted[j] += h;
            const auto rj = sys.residual(shifted);
            for (int i = 0; i < n; ++i) {
                const double fd = (rj[i] - r0[i]) / h;
                const double an = J(i, j);
                const double d = std::abs(fd - an) / std::max({std::abs(fd), std::abs(an), 1.0});
                if (d > worst) { worst = d; worst_i = i; worst_j = j; }
            }
        }
        BOOST_TEST_MESSAGE(sh.what << ": " << n << " unknowns, worst entry mismatch "
                           << 100 * worst << " % at (" << worst_i << "," << worst_j << ")");
        // A one-sided difference of phi carries its curvature; a wrong
        // derivative shows as tens of percent, not tenths of a permille.
        BOOST_CHECK_LT(worst, 1e-3);
    }
}

BOOST_AUTO_TEST_CASE(the_tree_on_a_network_by_three_routes)
{
    const auto deck_path = std::filesystem::path(__FILE__).parent_path() / "GROUPTREE.DATA";
    if (!std::filesystem::exists(deck_path)) {
        BOOST_TEST_MESSAGE("GROUPTREE.DATA not found, skipping");
        return;
    }
    const auto deck = Parser{}.parseFile(deck_path.string());
    const EclipseState es{deck};
    const Schedule schedule{deck, es};

    NetTreeCase tc;
    struct Case { const char* what; double target; std::array<double, 3> q0; bool fork; };
    const std::vector<Case> cases{
        {"target far above the tubing, chain",  2000.0, {400, 400, 300}, false},
        {"target holds every well, chain",       300.0, {400, 400, 300}, false},
        {"target between, chain",                600.0, {400, 400, 300}, false},
        {"target between, fork",                 600.0, {400, 400, 300}, true},
        {"one weak well, fork",                  600.0, {500, 500, 120}, true},
        {"one weak well, chain",                 600.0, {500, 500, 120}, false},
    };
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    for (const auto& c : cases) {
        for (const double g : {20.0, 30.0}) {
            const auto guess = NetTreeCase::guess(g);
            auto fb = tc.build(c.target, c.q0, c.fork);
            auto every = tc.build(c.target, c.q0, c.fork);
            every.setGroupActiveSet(true);
            auto outer = tc.build(c.target, c.q0, c.fork);
            const auto rf = NetworkSolve::solve(fb, guess, params, NetworkSolve::FullStep{});
            const auto re = NetworkSolve::solve(every, guess, params, NetworkSolve::FullStep{});
            const auto ro = NetworkSolve::solveWithTree(outer, guess, params, NetworkSolve::FullStep{});
            std::string sets;
            for (const auto& s : ro.sets) { sets += " " + s; }
            BOOST_TEST_MESSAGE(fmt::format(
                "{} (guess {} bar): outer {} passes / {} it{}{} sets{}; every-iteration {} it; fb {} it",
                c.what, g, ro.passes, ro.inner_iterations,
                ro.consistent ? "" : " NOT CONSISTENT", ro.cycled ? " CYCLED" : "", sets,
                re.converged ? std::to_string(re.iterations) : "FAILED",
                rf.converged ? std::to_string(rf.iterations) : "FAILED"));
            BOOST_TEST_MESSAGE("  outer" << rateList(outer, ro.result.well_rate)
                               << "  N1=" << convert::to(ro.result.node_pressure[1], bars)
                               << " N2=" << convert::to(ro.result.node_pressure[2], bars) << " bar");
            BOOST_CHECK(ro.result.converged);
            BOOST_CHECK(ro.consistent);
            BOOST_CHECK(!ro.cycled);
            BOOST_CHECK(re.converged);
            BOOST_CHECK(rf.converged);
            if (!ro.result.converged) { continue; }
            for (int w = 0; w < outer.numWells(); ++w) {
                if (re.converged) { BOOST_CHECK_CLOSE(ro.result.well_rate[w], re.well_rate[w], 0.5); }
                if (rf.converged) { BOOST_CHECK_CLOSE(ro.result.well_rate[w], rf.well_rate[w], 0.5); }
            }
            // The balancer, given the capacities at the converged pressures,
            // allocates what the Newton produced: its answer solves the network.
            const auto stein = steinAllocation(outer, ro.result.state, schedule);
            BOOST_REQUIRE(!stein.empty());
            std::string both;
            for (int w = 0; w < outer.numWells(); ++w) {
                both += fmt::format(" {}={:.1f}/{:.1f}", outer.wells()[w].name,
                                    ro.result.well_rate[w] * 86400.0, stein[w]);
                BOOST_CHECK_CLOSE(ro.result.well_rate[w] * 86400.0, stein[w], 0.5);
            }
            BOOST_TEST_MESSAGE("  newton/stein at the converged pressures:" << both);
        }
    }
}

// Constant phase fractions are exact when every phase line has the same
// shut-in pressure -- one cell, one pressure, fractions are mobility ratios
// -- and only then. With the water line's zero moved, the cut at the well's
// operating point and the cut where its own control puts it differ, and so
// do the two capacity measures. The outer loop is held to consistency under
// both; the answer is the Newton's either way, so what the fractions decide
// is only which set is tried.
BOOST_AUTO_TEST_CASE(constant_fractions_hold_only_on_a_common_shut_in)
{
    NetTreeCase tc;
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    using Sys = NetTreeCase::Sys;
    const auto guess = NetTreeCase::guess(20.0);
    for (const double water_shut : {120.0, 200.0, 80.0}) {
        std::array<double, 2> caps{};
        std::array<std::vector<double>, 2> rates;
        std::array<std::string, 2> sets;
        int i = 0;
        for (const auto frac : {Sys::CapacityFractions::Fixed, Sys::CapacityFractions::Ipr}) {
            auto sys = tc.build(700.0, {400, 400, 300}, true, Sys::Mode::Liquid, water_shut);
            sys.setCapacityFractions(frac);
            const auto r = NetworkSolve::solveWithTree(sys, guess, params, NetworkSolve::FullStep{});
            BOOST_CHECK(r.result.converged);
            BOOST_CHECK(r.consistent);
            // W1's capacity on the liquid mode at the converged pressures.
            sys.updateControls(r.result.state);
            const auto cw = sys.modeWeights(Sys::Mode::Liquid, {});
            caps[i] = sys.wellCapOnMode(r.result.state, 0, sys.ownAllowance(0), cw) * 86400.0;
            rates[i] = r.result.well_rate;
            sets[i] = r.sets.back();
            ++i;
        }
        BOOST_TEST_MESSAGE(fmt::format("water shut-in {} bar: W1 liquid capacity fixed {:.2f} / ipr {:.2f}"
                                       " sm3/d, sets {} / {}", water_shut, caps[0], caps[1], sets[0], sets[1]));
        if (water_shut == 120.0) {
            BOOST_CHECK_CLOSE(caps[0], caps[1], 1e-8);
        } else {
            BOOST_CHECK_GT(std::abs(caps[0] - caps[1]), 1e-3);
        }
        for (std::size_t w = 0; w < rates[0].size(); ++w) {
            BOOST_CHECK_CLOSE(rates[0][w], rates[1][w], 0.5);
        }
    }
}

// Looking for the flapping the outer loop could do: the balancer hands a set
// down, the network moves the pressures, the balancer takes it back. Fork and
// chain, targets from below every capacity to above their sum, stiff and
// soft wells, both starting guesses; count passes and any revisited set.
BOOST_AUTO_TEST_CASE(hunting_a_cycle_in_the_tree_pass)
{
    NetTreeCase tc;
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    int runs = 0, cycles = 0, inconsistent = 0, failed = 0, max_passes = 0;
    std::map<int, int> passes;
    long every_it = 0, outer_it = 0; int every_failed = 0;
    std::string worst;
    for (const bool fork : {false, true}) {
        for (const double scale : {0.5, 1.0, 2.0}) {
            for (double target = 100.0; target <= 1300.0; target += 50.0) {
                for (const double g : {20.0, 30.0}) {
                    const std::array<double, 3> q0{500.0 * scale, 400.0 * scale, 120.0 * scale};
                    auto sys = tc.build(target, q0, fork);
                    const auto r = NetworkSolve::solveWithTree(sys, NetTreeCase::guess(g), params,
                                                               NetworkSolve::FullStep{});
                    auto every = tc.build(target, q0, fork);
                    every.setGroupActiveSet(true);
                    const auto re = NetworkSolve::solve(every, NetTreeCase::guess(g), params,
                                                        NetworkSolve::FullStep{});
                    ++runs;
                    every_it += re.iterations; every_failed += re.converged ? 0 : 1;
                    outer_it += r.inner_iterations;
                    cycles += r.cycled ? 1 : 0;
                    inconsistent += (r.result.converged && !r.consistent) ? 1 : 0;
                    failed += r.result.converged ? 0 : 1;
                    ++passes[r.passes];
                    if (r.passes > max_passes || r.cycled) {
                        max_passes = std::max(max_passes, r.passes);
                        std::string sets;
                        for (const auto& s : r.sets) { sets += " " + s; }
                        worst = fmt::format("{} target {} scale {} guess {}: {} passes{}{}",
                                            fork ? "fork" : "chain", target, scale, g, r.passes,
                                            r.cycled ? " CYCLED" : "", sets);
                    }
                }
            }
        }
    }
    std::string hist;
    for (const auto& [n, k] : passes) { hist += fmt::format(" {}:{}", n, k); }
    BOOST_TEST_MESSAGE(fmt::format("{} runs: {} cycled, {} inconsistent, {} failed; passes{}; "
                                   "newton it outer {} vs every-iteration {} ({} failed)",
                                   runs, cycles, inconsistent, failed, hist, outer_it, every_it, every_failed));
    BOOST_TEST_MESSAGE("longest: " << worst);
    BOOST_CHECK_EQUAL(cycles, 0);
    BOOST_CHECK_EQUAL(inconsistent, 0);
    BOOST_CHECK_EQUAL(failed, 0);
}

// The reduced form: Newton on the node pressures alone, the tree walk as the
// well model (NetworkReducedSolve.hpp), against the outer loop on the same
// cases and the same sweep. Same fixed points by construction; what differs
// is how the kinks are met -- crossed inside a differenced Jacobian and a
// line search here, frozen out of the Newton there -- and what each costs.
BOOST_AUTO_TEST_CASE(the_reduced_form_against_the_outer_loop)
{
    NetTreeCase tc;
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    struct Case { const char* what; double target; std::array<double, 3> q0; bool fork; };
    const std::vector<Case> cases{
        {"target far above the tubing, chain",  2000.0, {400, 400, 300}, false},
        {"target holds every well, chain",       300.0, {400, 400, 300}, false},
        {"target between, fork",                 600.0, {400, 400, 300}, true},
        {"one weak well, fork",                  600.0, {500, 500, 120}, true},
        {"one weak well, chain",                 600.0, {500, 500, 120}, false},
    };
    for (const auto& c : cases) {
        for (const double g : {20.0, 30.0}) {
            const auto guess = NetTreeCase::guess(g);
            auto outer = tc.build(c.target, c.q0, c.fork);
            auto reduced = tc.build(c.target, c.q0, c.fork);
            outer.resetLookups(); reduced.resetLookups();
            const auto ro = NetworkSolve::solveWithTree(outer, guess, params, NetworkSolve::FullStep{});
            const auto rr = NetworkSolve::solveReduced(reduced, guess, params);
            BOOST_TEST_MESSAGE(fmt::format(
                "{} (guess {} bar): reduced {} in {} it, {} evaluations, {} stalls, {} set changes,"
                " {} lookups [{}]; outer {} passes / {} it, {} lookups",
                c.what, g, rr.converged ? "converged" : "FAILED", rr.iterations, rr.evaluations,
                rr.stalls, rr.set_changes, reduced.lookups(), rr.sets,
                ro.passes, ro.inner_iterations, outer.lookups()));
            BOOST_TEST_MESSAGE("  reduced" << rateList(reduced, rr.well_rate)
                               << "  N1=" << convert::to(rr.node_pressure[1], bars)
                               << " N2=" << convert::to(rr.node_pressure[2], bars) << " bar");
            BOOST_CHECK(rr.converged);
            BOOST_CHECK(ro.result.converged);
            if (!rr.converged || !ro.result.converged) { continue; }
            for (int w = 0; w < outer.numWells(); ++w) {
                BOOST_CHECK_CLOSE(rr.well_rate[w], ro.result.well_rate[w], 0.5);
            }
            for (int nd = 1; nd <= outer.numNodes(); ++nd) {
                BOOST_CHECK_CLOSE(rr.node_pressure[nd], ro.result.node_pressure[nd], 0.5);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(the_reduced_form_on_the_sweep)
{
    NetTreeCase tc;
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    int runs = 0, failed = 0, stalls = 0, disagree = 0, max_it = 0;
    long it_reduced = 0, it_outer = 0, eval = 0, look_reduced = 0, look_outer = 0, set_changes = 0;
    std::map<int, int> hist;
    std::string worst;
    for (const bool fork : {false, true}) {
        for (const double scale : {0.5, 1.0, 2.0}) {
            for (double target = 100.0; target <= 1300.0; target += 50.0) {
                for (const double g : {20.0, 30.0}) {
                    const std::array<double, 3> q0{500.0 * scale, 400.0 * scale, 120.0 * scale};
                    auto outer = tc.build(target, q0, fork);
                    auto reduced = tc.build(target, q0, fork);
                    outer.resetLookups(); reduced.resetLookups();
                    const auto ro = NetworkSolve::solveWithTree(outer, NetTreeCase::guess(g), params,
                                                                NetworkSolve::FullStep{});
                    const auto rr = NetworkSolve::solveReduced(reduced, NetTreeCase::guess(g), params);
                    ++runs;
                    it_reduced += rr.iterations; it_outer += ro.inner_iterations;
                    eval += rr.evaluations; stalls += rr.stalls; set_changes += rr.set_changes;
                    look_reduced += reduced.lookups(); look_outer += outer.lookups();
                    ++hist[rr.iterations];
                    failed += rr.converged ? 0 : 1;
                    bool off = false;
                    if (rr.converged && ro.result.converged) {
                        for (int w = 0; w < outer.numWells(); ++w) {
                            const double a = rr.well_rate[w], b = ro.result.well_rate[w];
                            if (std::abs(a - b) > 0.005 * std::max(std::abs(b), 1.0 / 86400.0)) { off = true; }
                        }
                    }
                    disagree += off ? 1 : 0;
                    if (rr.iterations > max_it || !rr.converged || off) {
                        max_it = std::max(max_it, rr.iterations);
                        worst = fmt::format("{} target {} scale {} guess {}: {} it, {} stalls{}{} [{}]"
                                            " reduced{} outer{}",
                                            fork ? "fork" : "chain", target, scale, g, rr.iterations,
                                            rr.stalls, rr.converged ? "" : " FAILED",
                                            off ? " DISAGREES" : "", rr.sets,
                                            rateList(reduced, rr.well_rate),
                                            rateList(outer, ro.result.well_rate));
                    }
                }
            }
        }
    }
    std::string h;
    for (const auto& [k, v] : hist) { h += fmt::format(" {}:{}", k, v); }
    BOOST_TEST_MESSAGE(fmt::format(
        "{} runs: reduced {} failed, {} disagree, {} stalls, {} set changes, {} it, {} evaluations,"
        " {} lookups; outer {} it, {} lookups; reduced iterations{}",
        runs, failed, disagree, stalls, set_changes, it_reduced, eval, look_reduced, it_outer,
        look_outer, h));
    BOOST_TEST_MESSAGE("worst: " << worst);
    BOOST_CHECK_EQUAL(failed, 0);
    BOOST_CHECK_EQUAL(disagree, 0);
}

// The reduced Jacobian two ways: differenced along p, and eliminated from
// the full system's assembled one. They must agree wherever the differences
// stay on one piece -- and then the eliminated one is the derivative, exact,
// for one assembly and N small solves instead of N + 1 evaluations.
BOOST_AUTO_TEST_CASE(the_reduced_jacobian_by_elimination)
{
    NetTreeCase tc;
    struct Case { const char* what; double target; std::array<double, 3> q0; bool fork; };
    const std::vector<Case> cases{
        {"target far above the tubing, chain",  2000.0, {400, 400, 300}, false},
        {"target holds every well, chain",       300.0, {400, 400, 300}, false},
        {"one weak well, fork",                  600.0, {500, 500, 120}, true},
        {"one weak well, chain",                 600.0, {500, 500, 120}, false},
    };
    for (const auto& c : cases) {
        auto sys = tc.build(c.target, c.q0, c.fork);
        sys.setGroupActiveSet(true);
        sys.setExactPotential(true);
        // Where the wells are on their tubing, so the node rows feel the wells:
        // at the guess pressures they are on bhp or share and the Jacobian is
        // the trivial one.
        const std::vector<double> p{convert::from(20.0, bars), convert::from(34.0, bars),
                                    convert::from(41.0, bars)};
        const auto r = sys.reducedResidual(p);
        const auto S = NetworkSolve::reducedJacobianByElimination(sys);
        const int n = sys.numNodes();
        double worst = 0.0;
        std::string entries;
        const double h = 1e-3 * unit::barsa;
        for (int j = 0; j < n; ++j) {
            auto pj = p; pj[j + 1] += h;
            const auto rj = sys.reducedResidual(pj);
            for (int i = 0; i < n; ++i) {
                const double fd = (rj[i] - r[i]) / h * unit::barsa;   // per bar
                const double an = S(i, j) * unit::barsa;
                entries += fmt::format(" [{},{}] {:.4f}/{:.4f}", i, j, fd, an);
                worst = std::max(worst, std::abs(fd - an) / std::max({std::abs(fd), std::abs(an), 1e-2}));
            }
        }
        sys.reducedResidual(p);       // leave the state where S was taken
        BOOST_TEST_MESSAGE(c.what << " [" << sys.treeSignature() << "]: fd/eliminated" << entries
                           << fmt::format("  worst {:.3g} %", 100 * worst));
        BOOST_CHECK_LT(worst, 1e-3);
    }

    // And the sweep with the eliminated Jacobian, against the differenced.
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    int runs = 0, failed = 0, stalls = 0, disagree = 0;
    long it_e = 0, it_d = 0, eval_e = 0, eval_d = 0, look_e = 0, look_d = 0;
    for (const bool fork : {false, true}) {
        for (const double scale : {0.5, 1.0, 2.0}) {
            for (double target = 100.0; target <= 1300.0; target += 50.0) {
                for (const double g : {20.0, 30.0}) {
                    const std::array<double, 3> q0{500.0 * scale, 400.0 * scale, 120.0 * scale};
                    auto a = tc.build(target, q0, fork);
                    auto b = tc.build(target, q0, fork);
                    a.resetLookups(); b.resetLookups();
                    const auto re = NetworkSolve::solveReduced(a, NetTreeCase::guess(g), params, true);
                    const auto rd = NetworkSolve::solveReduced(b, NetTreeCase::guess(g), params, false);
                    ++runs;
                    failed += re.converged ? 0 : 1;
                    stalls += re.stalls;
                    it_e += re.iterations; it_d += rd.iterations;
                    eval_e += re.evaluations; eval_d += rd.evaluations;
                    look_e += a.lookups(); look_d += b.lookups();
                    if (re.converged && rd.converged) {
                        for (int w = 0; w < a.numWells(); ++w) {
                            if (std::abs(re.well_rate[w] - rd.well_rate[w])
                                > 0.005 * std::max(std::abs(rd.well_rate[w]), 1.0 / 86400.0)) {
                                ++disagree; break;
                            }
                        }
                    }
                }
            }
        }
    }
    BOOST_TEST_MESSAGE(fmt::format(
        "{} runs, eliminated Jacobian: {} failed, {} stalls, {} disagree; it {} vs {} differenced;"
        " evaluations {} vs {}; lookups {} vs {}",
        runs, failed, stalls, disagree, it_e, it_d, eval_e, eval_d, look_e, look_d));
    BOOST_CHECK_EQUAL(failed, 0);
    BOOST_CHECK_EQUAL(disagree, 0);
}

BOOST_AUTO_TEST_CASE(model5_both_trees_from_the_deck)
{
    const std::string path = kNetworkDecks + "NETWORK_MODEL5_STDW_AUTOCHK.DATA";
    if (!std::filesystem::exists(path)) {
        BOOST_TEST_MESSAGE("opm-tests not present, skipping");
        return;
    }
    DeckTrees dt(path);
    const int step = 1;                         // after WELOPEN: B-1H, B-2H, B-3H, C-1H
    auto sys = dt.build(step, {});
    BOOST_TEST_MESSAGE("MODEL5:" << dt.description);
    // The two trees differ: C1's group parent M5N is not a node.
    BOOST_CHECK(dt.gidx.count("M5N"));
    BOOST_CHECK(std::find(dt.node_order.begin(), dt.node_order.end(), "M5N") == dt.node_order.end());
    BOOST_CHECK_EQUAL(sys.numWells(), 4);
    for (const double g : {21.0, 30.0}) {
        threeRoutesAndOracle(dt, {}, step, g, fmt::format("model5, deck as is, guess {}", g).c_str());
    }
    // Stiffer and softer wells, and a tighter B1.
    for (const double j : {1.0, 4.0}) {
        DeckTrees::Ipr ipr; ipr.j_scale = j;
        threeRoutesAndOracle(dt, ipr, step, 21.0, fmt::format("model5, j_scale {}", j).c_str());
    }
}

BOOST_AUTO_TEST_CASE(network01_both_trees_from_the_deck)
{
    const std::string path = kNetworkDecks + "NETWORK-01.DATA";
    if (!std::filesystem::exists(path)) {
        BOOST_TEST_MESSAGE("opm-tests not present, skipping");
        return;
    }
    DeckTrees dt(path);
    const int step = 1;
    auto sys = dt.build(step, {});
    BOOST_TEST_MESSAGE("NETWORK-01:" << dt.description);
    BOOST_CHECK_EQUAL(sys.numWells(), 2);
    for (const double g : {80.0, 90.0}) {
        threeRoutesAndOracle(dt, {}, step, g, fmt::format("network-01, deck as is, guess {}", g).c_str());
    }
    DeckTrees::Ipr ipr; ipr.j_scale = 4.0;
    threeRoutesAndOracle(dt, ipr, step, 80.0, "network-01, j_scale 4");
}

// The same sweep as on the hand-built network, on MODEL5: the B1 target
// from below what one well makes to above what three make, three well
// stiffnesses, two starting pressures.
BOOST_AUTO_TEST_CASE(model5_sweep)
{
    const std::string path = kNetworkDecks + "NETWORK_MODEL5_STDW_AUTOCHK.DATA";
    if (!std::filesystem::exists(path)) {
        BOOST_TEST_MESSAGE("opm-tests not present, skipping");
        return;
    }
    DeckTrees dt(path);
    const int step = 1;
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    int runs = 0, failed = 0, cycled = 0, inconsistent = 0, disagree = 0, r_failed = 0, r_stalls = 0;
    long it_outer = 0, it_reduced = 0;
    std::map<int, int> passes;
    std::string worst;
    for (const double j : {1.0, 2.0, 4.0}) {
        for (double target = 1000.0; target <= 9000.0; target += 500.0) {
            for (const double g : {21.0, 30.0}) {
                DeckTrees::Ipr ipr; ipr.j_scale = j;
                auto outer = dt.build(step, ipr);
                auto reduced = dt.build(step, ipr);
                const int b1 = dt.gidx.at("B1");
                outer.setGroupLimit(b1, DeckTrees::Sys::Mode::Oil, target / 86400.0);
                reduced.setGroupLimit(b1, DeckTrees::Sys::Mode::Oil, target / 86400.0);
                const auto ro = NetworkSolve::solveWithTree(outer, dt.guess(g), params, NetworkSolve::FullStep{});
                const auto rr = NetworkSolve::solveReduced(reduced, dt.guess(g), params, true);
                ++runs;
                failed += ro.result.converged ? 0 : 1;
                cycled += ro.cycled ? 1 : 0;
                inconsistent += (ro.result.converged && !ro.consistent) ? 1 : 0;
                r_failed += rr.converged ? 0 : 1;
                r_stalls += rr.stalls;
                it_outer += ro.inner_iterations; it_reduced += rr.iterations;
                ++passes[ro.passes];
                bool off = false;
                if (ro.result.converged && rr.converged) {
                    for (int w = 0; w < outer.numWells(); ++w) {
                        if (std::abs(ro.result.well_rate[w] - rr.well_rate[w])
                            > 0.005 * std::max(std::abs(rr.well_rate[w]), 1.0 / 86400.0)) { off = true; }
                    }
                }
                disagree += off ? 1 : 0;
                if (!ro.result.converged || ro.cycled || off || !rr.converged || rr.stalls) {
                    std::string sets;
                    for (const auto& s : ro.sets) { sets += " " + s; }
                    worst = fmt::format("j {} target {} guess {}: outer {} passes{}{}{}, reduced {} it {} stalls{}"
                                        " [{}] outer{} reduced{}", j, target, g, ro.passes,
                                        ro.result.converged ? "" : " FAILED", ro.cycled ? " CYCLED" : "",
                                        off ? " DISAGREES" : "", rr.iterations, rr.stalls,
                                        rr.converged ? "" : " FAILED", sets,
                                        rateList(outer, ro.result.well_rate), rateList(reduced, rr.well_rate));
                }
            }
        }
    }
    std::string hist;
    for (const auto& [n, k] : passes) { hist += fmt::format(" {}:{}", n, k); }
    BOOST_TEST_MESSAGE(fmt::format("MODEL5, {} runs: outer {} failed, {} cycled, {} inconsistent, passes{}, {} it;"
                                   " reduced {} failed, {} stalls, {} it; {} disagree",
                                   runs, failed, cycled, inconsistent, hist, it_outer,
                                   r_failed, r_stalls, it_reduced, disagree));
    if (!worst.empty()) { BOOST_TEST_MESSAGE("worst: " << worst); }
    BOOST_CHECK_EQUAL(failed, 0);
    BOOST_CHECK_EQUAL(cycled, 0);
    BOOST_CHECK_EQUAL(inconsistent, 0);
    BOOST_CHECK_EQUAL(disagree, 0);
}

BOOST_AUTO_TEST_CASE(generated_two_tree_cases)
{
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    struct Size { int wells, nodes, groups, seeds; };
    const std::vector<Size> sizes{{10, 4, 4, 6}, {50, 15, 12, 4}, {200, 40, 30, 2}};
    for (const auto& sz : sizes) {
        int runs = 0, failed = 0, cycled = 0, inconsistent = 0, r_failed = 0, r_stalls = 0,
            disagree = 0, two_solutions = 0, oracle_off = 0, oracle_violates = 0, oracle_rejected = 0,
            max_passes = 0, off_branch_total = 0, off_axis = 0;
        long it_outer = 0, it_reduced = 0, look_outer = 0, look_reduced = 0;
        double ms_outer = 0.0, ms_reduced = 0.0, ms_oracle = 0.0;
        for (int seed = 1; seed <= sz.seeds; ++seed) {
            std::string what;
            const auto text = generateDeck({sz.wells, sz.nodes, sz.groups, static_cast<unsigned>(seed)}, &what);
            DeckTrees dt(text, DeckTrees::FromText{});
            for (const double j : {1.0, 2.5}) {
                DeckTrees::Ipr ipr; ipr.j_scale = j;
                auto outer = dt.build(0, ipr);
                auto reduced = dt.build(0, ipr);
                outer.resetLookups(); reduced.resetLookups();
                const auto guess = dt.guess(20.0);
                using clock = std::chrono::steady_clock;
                auto t0 = clock::now();
                const auto ro = NetworkSolve::solveWithTree(outer, guess, params, NetworkSolve::FullStep{});
                auto t1 = clock::now();
                const auto rr = NetworkSolve::solveReduced(reduced, guess, params, true);
                auto t2 = clock::now();
                ms_outer += std::chrono::duration<double, std::milli>(t1 - t0).count();
                ms_reduced += std::chrono::duration<double, std::milli>(t2 - t1).count();
                ++runs;
                failed += ro.result.converged ? 0 : 1;
                cycled += ro.cycled ? 1 : 0;
                inconsistent += (ro.result.converged && !ro.consistent) ? 1 : 0;
                r_failed += rr.converged ? 0 : 1;
                r_stalls += rr.stalls;
                off_axis += (ro.result.off_axis > 0 ? 1 : 0) + (rr.off_axis > 0 ? 1 : 0);
                it_outer += ro.inner_iterations; it_reduced += rr.iterations;
                look_outer += outer.lookups(); look_reduced += reduced.lookups();
                max_passes = std::max(max_passes, ro.passes);
                bool off = false;
                if (ro.result.converged && rr.converged) {
                    for (int w = 0; w < outer.numWells(); ++w) {
                        if (std::abs(ro.result.well_rate[w] - rr.well_rate[w])
                            > 0.005 * std::max(std::abs(rr.well_rate[w]), 1.0 / 86400.0)) { off = true; }
                    }
                }
                // Two different answers, each self-consistent by the walk, are
                // two fixed points, not a fault of either route.
                if (off && ro.consistent) { ++two_solutions; off = false; }
                disagree += off ? 1 : 0;
                // The oracle against the reduced answer, whose wells are on the
                // stable crossing by construction; and how many of the outer
                // loop's thp wells are not -- its Newton can land on the hump.
                double worst_gap = 0.0, worst_gap_outer = 0.0;
                bool rejected = false;
                int off_branch = 0;
                if (rr.converged) {
                    auto t3 = clock::now();
                    const auto stein = steinAllocation(reduced, reduced.reducedState(), *dt.schedule, 0);
                    ms_oracle += std::chrono::duration<double, std::milli>(clock::now() - t3).count();
                    if (stein.empty()) { rejected = true; ++oracle_rejected; }
                    else {
                        for (int w = 0; w < reduced.numWells(); ++w) {
                            const double mine = rr.well_rate[w] * 86400.0;
                            worst_gap = std::max(worst_gap, std::abs(mine - stein[w]) / std::max(std::abs(stein[w]), 1.0));
                        }
                        // An oracle allocation above some group's own limit is
                        // the balancer walking through a transparent group
                        // (a_nested_target_under_a_transparent_group), not a
                        // disagreement about the answer.
                        bool violates = false;
                        const auto& groups = reduced.groups();
                        for (int g = 0; g < reduced.numGroups() && !violates; ++g) {
                            if (!(groups[g].target > 0.0)) { continue; }
                            const auto c = DeckTrees::Sys::modeWeights(groups[g].mode, {});
                            double on = 0.0;
                            for (int w = 0; w < reduced.numWells(); ++w) {
                                bool under = false;
                                for (int a = reduced.wells()[w].group; a >= 0; a = groups[a].parent) { if (a == g) { under = true; break; } }
                                if (!under) { continue; }
                                const auto& well = reduced.wells()[w];
                                const double q_oil = stein[w] / 86400.0;
                                const double bhp = (q_oil - well.ipr_a[1]) / well.ipr_b[1];
                                for (int ph = 0; ph < DeckTrees::Sys::NP; ++ph) { on += c[ph] * std::max(well.ipr_a[ph] + well.ipr_b[ph] * bhp, 0.0); }
                            }
                            if (on > groups[g].target * 1.005) { violates = true; }
                        }
                        if (worst_gap > 0.005) { (violates ? oracle_violates : oracle_off) += 1; }
                    }
                }
                if (ro.result.converged) {
                    const auto stein = steinAllocation(outer, ro.result.state, *dt.schedule, 0);
                    for (int w = 0; w < outer.numWells() && !stein.empty(); ++w) {
                        const double mine = ro.result.well_rate[w] * 86400.0;
                        worst_gap_outer = std::max(worst_gap_outer, std::abs(mine - stein[w]) / std::max(std::abs(stein[w]), 1.0));
                    }
                    for (int w = 0; w < outer.numWells(); ++w) {
                        if (outer.control(w) != DeckTrees::Sys::Control::Thp) { continue; }
                        const auto& well = outer.wells()[w];
                        const double p_node = ro.result.node_pressure[well.node];
                        const double q_cross = outer.thpPotential(well, p_node);
                        if (q_cross > 0.0 && q_cross < std::numeric_limits<double>::max()
                            && std::abs(q_cross - ro.result.well_rate[w]) > 0.01 * q_cross) { ++off_branch; }
                    }
                    off_branch_total += off_branch;
                }
                if (ro.result.off_axis || rr.off_axis) {
                    BOOST_TEST_MESSAGE(fmt::format("  seed {} j {}: off the tables -- outer {} ({}), reduced {} ({})",
                                                   seed, j, ro.result.off_axis, outer.offAxisNote(),
                                                   rr.off_axis, reduced.offAxisNote()));
                }
                if (!ro.result.converged || ro.cycled || !ro.consistent || off || !rr.converged
                    || rr.stalls || rejected || worst_gap > 0.005 || off_branch) {
                    std::string sets;
                    for (const auto& s : ro.sets) { sets += " " + s; }
                    BOOST_TEST_MESSAGE(fmt::format(
                        "  seed {} j {}: {}; outer {} passes / {} it{}{}{}, {} thp wells off the stable crossing,"
                        " oracle gap {:.3g} %; reduced {} it {} stalls{}, oracle {}{}",
                        seed, j, what, ro.passes, ro.inner_iterations,
                        ro.result.converged ? "" : " FAILED", ro.cycled ? " CYCLED" : "",
                        ro.consistent ? "" : " INCONSISTENT", off_branch, 100 * worst_gap_outer,
                        rr.iterations, rr.stalls, rr.converged ? "" : " FAILED",
                        rejected ? "rejected" : fmt::format("gap {:.3g} %", 100 * worst_gap),
                        off ? "; outer/reduced DISAGREE" : "")
                        + (sz.wells <= 10 ? " sets" + sets : std::string{}));
                }
            }
        }
        BOOST_TEST_MESSAGE(fmt::format(
            "{}/{}/{} wells/nodes/groups, {} runs: outer {} failed, {} cycled, {} inconsistent, max {} passes,"
            " {} it, {} lookups, {:.1f} ms/run; reduced {} failed, {} stalls, {} it, {} lookups, {:.1f} ms/run;"
            " {} disagree, {} two self-consistent answers, {} answers off the tables; {} outer thp wells off"
            " the stable crossing; oracle vs reduced {} rejected, {} off, {} above a group's own limit, {:.1f} ms/run",
            sz.wells, sz.nodes, sz.groups, runs, failed, cycled, inconsistent, max_passes, it_outer,
            look_outer, ms_outer / runs, r_failed, r_stalls, it_reduced, look_reduced, ms_reduced / runs,
            disagree, two_solutions, off_axis, off_branch_total, oracle_rejected, oracle_off, oracle_violates,
            ms_oracle / std::max(runs, 1)));
        BOOST_CHECK_EQUAL(failed, 0);
        BOOST_CHECK_EQUAL(cycled, 0);
        BOOST_CHECK_EQUAL(inconsistent, 0);
        BOOST_CHECK_EQUAL(r_failed, 0);
        BOOST_CHECK_EQUAL(disagree, 0);
        BOOST_CHECK_EQUAL(oracle_rejected, 0);
        BOOST_CHECK_EQUAL(oracle_off, 0);
        BOOST_CHECK_EQUAL(off_axis, 0);
    }
}

// The guard itself: the instance whose tables were sized for ten wells, at
// two hundred. Both routes converge, and both say so -- the answer needed
// lookups off the tables' axes, so it is an answer to nothing.
BOOST_AUTO_TEST_CASE(an_answer_off_the_tables_says_so)
{
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    std::string what;
    const auto text = generateDeck({200, 40, 30, 2u, /*size_tables=*/false}, &what);
    DeckTrees dt(text, DeckTrees::FromText{});
    DeckTrees::Ipr ipr; ipr.j_scale = 2.5;
    auto reduced = dt.build(0, ipr);
    auto outer = dt.build(0, ipr);
    const auto rr = NetworkSolve::solveReduced(reduced, dt.guess(20.0), params, true);
    const auto ro = NetworkSolve::solveWithTree(outer, dt.guess(20.0), params, NetworkSolve::FullStep{});
    BOOST_TEST_MESSAGE(fmt::format("{}: reduced {} with {} lookups off the tables ({}); outer {} with {} ({})",
                                   what, rr.converged ? "converged" : "FAILED", rr.off_axis, reduced.offAxisNote(),
                                   ro.result.converged ? "converged" : "FAILED", ro.result.off_axis, outer.offAxisNote()));
    BOOST_CHECK(rr.converged);
    BOOST_CHECK_GT(rr.off_axis, 0);
    BOOST_CHECK(ro.result.converged);
    BOOST_CHECK_GT(ro.result.off_axis, 0);
}

// The smallest generated instance that goes wrong, in full: the two trees,
// what each well may make, what each route and the oracle gave it, and --
// for a singular first Jacobian -- which bound group has no held child.
// A group's own target two levels under a bound parent, with a transparent
// group between: seed 3 of the generator in three wells. PLAT wants 3000;
// G5 may make 500 of it, W1 2000, so PLAT cannot be met and G5 must stay at
// 500. The oracle is held to it as well.
BOOST_AUTO_TEST_CASE(a_nested_target_under_a_transparent_group)
{
    std::string d = deckPreamble(3, 5, 2);
    d += "BRANPROP\n 'G1' 'PLAT' 9999 1* /\n/\nNODEPROP\n 'PLAT' 20.0 NO NO 1* /\n 'G1' 1* NO NO 1* /\n/\n";
    d += "GRUPTREE\n 'PLAT' 'FIELD' /\n 'G1' 'PLAT' /\n 'G3' 'PLAT' /\n 'G5' 'G3' /\n/\n";
    d += "WELSPECS\n 'W1' 'G1' 1 1 7000 'OIL' /\n 'W2' 'G5' 2 1 7000 'OIL' /\n 'W3' 'G5' 3 1 7000 'OIL' /\n/\n";
    d += "COMPDAT\n 'W1' 1 1 1 2 'OPEN' 1* 1* 0.2 /\n 'W2' 2 1 1 2 'OPEN' 1* 1* 0.2 /\n 'W3' 3 1 1 2 'OPEN' 1* 1* 0.2 /\n/\n";
    d += "WCONPROD\n 'W1' 'OPEN' 'GRUP' 2000 4* 100 /\n 'W2' 'OPEN' 'GRUP' 2000 4* 100 /\n 'W3' 'OPEN' 'GRUP' 2000 4* 100 /\n/\n";
    d += "GCONPROD\n 'PLAT' 'ORAT' 3000 3* 'RATE' /\n 'G5' 'ORAT' 500 3* 'RATE' /\n/\nTSTEP\n1 /\nEND\n";
    DeckTrees dt(d, DeckTrees::FromText{});
    auto sys = dt.build(0, {});
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    const auto rr = NetworkSolve::solveReduced(sys, dt.guess(20.0), params, true);
    BOOST_REQUIRE(rr.converged);
    const auto stein = steinAllocation(sys, sys.reducedState(), *dt.schedule, 0);
    BOOST_REQUIRE(!stein.empty());
    std::string both;
    double g5_mine = 0.0, g5_stein = 0.0;
    for (int w = 0; w < sys.numWells(); ++w) {
        both += fmt::format(" {}={:.1f}/{:.1f}", sys.wells()[w].name, rr.well_rate[w] * 86400.0, stein[w]);
        if (w > 0) { g5_mine += rr.well_rate[w] * 86400.0; g5_stein += stein[w]; }
    }
    BOOST_TEST_MESSAGE("reduced/stein:" << both << fmt::format("  G5 total {:.1f} / {:.1f} against its 500",
                                                                g5_mine, g5_stein));
    BOOST_CHECK_CLOSE(g5_mine, 500.0, 0.5);
    BOOST_CHECK_CLOSE(rr.well_rate[0] * 86400.0, 2000.0, 0.5);
    // The balancer, as vendored, hands G5 2000: a group's own limit under a
    // group without a guide rate is walked through. Recorded, not enforced.
    BOOST_WARN_LE(g5_stein, 500.0 * 1.005);
}

BOOST_AUTO_TEST_SUITE_END()
