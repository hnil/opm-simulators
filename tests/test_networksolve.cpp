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

// The production system: its rows, its controls, the analytic Jacobian,
// and the dumped simulator failures they were written from.

BOOST_AUTO_TEST_SUITE(NetworkSolveBench)

// Prototype: the same formulation on a production network.
//
// A rate becomes three numbers instead of one, and that is the whole difference.
// VFPPROD works the water and gas fractions out of the triple it is given, so
// they never become unknowns, and mixing at a node is a sum. The Newton, the
// active set and the scaling are the injection ones, unchanged.
//
// PROD -> FIELD, terminal at 80 bar, two producers on one node.
BOOST_AUTO_TEST_CASE(production_network_prototype)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    const auto sm3d = cubic(meter) / day;

    const auto deck = Parser{}.parseString(vfp_prod);
    const VFPProdTable table(deck["VFPPROD"].front(), /*gaslift_opt_active=*/false, UnitSystem{});
    VFPProdProperties<double> props;
    props.addTable(table);
    const UnitSystem units{};

    Sys system(props, units);
    system.setTerminalPressure(convert::from(80.0, bars));
    system.addNode(NetworkSolve::Node{"FIELD", -1, NetworkSolve::NoTable}, 0.0);
    system.addNode(NetworkSolve::Node{"PROD", 0, 3}, 0.0);

    // Two producers, water-cut about 0.3, GOR near the table's single value.
    for (const auto& [name, productivity] : std::initializer_list<std::pair<const char*, double>>{
             {"P-1", 1.0}, {"P-2", 0.7}}) {
        Sys::Well w;
        w.name = name;
        w.node = 1;
        w.vfp_table = 3;
        w.bhp_limit = convert::from(40.0, bars);
        w.oil_rate_limit = convert::from(600.0, sm3d);
        // q_p = a_p - b_p * bhp: production falls as bhp rises.
        const double q0 = convert::from(400.0 * productivity, sm3d);
        const double slope = q0 / convert::from(120.0, bars);
        for (int ph = 0; ph < Sys::NP; ++ph) {
            const double share = (ph == 0) ? 0.3 : (ph == 1) ? 0.7 : 70.0;  // water, oil, gas
            w.ipr_a[ph] = share * q0 * 2.0;
            w.ipr_b[ph] = -share * slope * 2.0;
        }
        system.addWell(w);
    }
    system.finish();

    BOOST_TEST_MESSAGE("unknowns: " << system.size() << "  (nodes " << system.numNodes()
                       << ", wells " << system.numWells() << ")");

    const std::vector<double> guess{convert::from(80.0, bars), convert::from(90.0, bars)};
    const auto r = NetworkSolve::solve(system, guess, kParams, NetworkSolve::FullStep{});
    BOOST_TEST_MESSAGE((r.converged ? "converged in " : "FAILED after ") << r.iterations
                       << " iterations, residual " << r.residual
                       << (r.control_trace.empty() ? "" : "  controls " + r.control_trace));
    BOOST_REQUIRE(r.converged);

    const double p_prod = r.node_pressure[1];
    BOOST_TEST_MESSAGE("PROD node " << convert::to(p_prod, bars) << " bar, oil "
                       << convert::to(r.well_rate[0], sm3d) << " + "
                       << convert::to(r.well_rate[1], sm3d) << " sm3/d");

    // The node sits above its terminal, and every well is producing.
    BOOST_CHECK_GT(convert::to(p_prod, bars), 80.0);
    for (const double q : r.well_rate) {
        BOOST_CHECK_GT(convert::to(q, sm3d), 0.0);
    }
}

// The group equations against the rule they replace.
//
// Given linearised wells and a set of node pressures, the old way to place a
// group's rate is arithmetic: cap each well by what it can actually take, share
// the target out by guide rate, and whenever a well's share exceeds its cap, fix
// it there, take it out of the pool and share the remainder among the rest. The
// new way is two equations and a multiplier. They must agree, and if they do not
// the equations are wrong -- this is the check that they are not.
BOOST_AUTO_TEST_CASE(group_equations_match_the_rule_based_allocation)
{
    const auto sm3d = cubic(meter) / day;

    // The rule, worked out at whatever pressures the equations settled on: cap
    // each well by what it can actually take, share the target by guide rate,
    // and whenever a share exceeds a cap, fix that well there, drop it from the
    // pool and share the remainder among the rest.
    auto ruleBased = [&](const NetworkSolve::InjectionSystem<double>& system,
                         const std::vector<double>& node_pressure,
                         const double target) {
        const auto& wells = system.wells();
        const int n = static_cast<int>(wells.size());
        std::vector<double> cap(n), share(n, 0.0);
        std::vector<bool> pooled(n, true);
        for (int w = 0; w < n; ++w) {
            const double p = node_pressure[wells[w].node];
            cap[w] = std::min({system.thpPotential(wells[w], p, /*cap_by_rate_limit=*/true), wells[w].rate_limit,
                               NetworkSolve::InjectionSystem<double>::ipr(wells[w], wells[w].bhp_limit)});
        }
        double remaining = target;
        for (int pass = 0; pass <= n; ++pass) {
            double guides = 0.0;
            for (int w = 0; w < n; ++w) {
                if (pooled[w]) {
                    guides += wells[w].guide;
                }
            }
            if (guides <= 0.0) {
                break;
            }
            bool fixed_one = false;
            for (int w = 0; w < n; ++w) {
                if (!pooled[w]) {
                    continue;
                }
                const double s = remaining * wells[w].guide / guides;
                if (s > cap[w]) {
                    share[w] = cap[w];
                    pooled[w] = false;
                    remaining -= cap[w];
                    fixed_one = true;
                    break;
                }
                share[w] = s;
            }
            if (!fixed_one) {
                break;
            }
        }
        return share;
    };

    auto compare = [&](const char* what, NetworkCase& c, const double target,
                       const bool required) {
        auto system = c.system();
        const auto r = NetworkSolve::solve(system, c.nodePressures(kStart), kParams, NetworkSolve::FullStep{});
        if (!r.converged) {
            BOOST_TEST_MESSAGE(what << ": the solve does not converge, so the equations cannot "
                               "be compared here yet");
            BOOST_CHECK(!required);
            return;
        }
        const auto share = ruleBased(system, r.node_pressure, target);
        double rule_total = 0.0, solved_total = 0.0;
        for (std::size_t w = 0; w < share.size(); ++w) {
            BOOST_TEST_MESSAGE("  " << system.wells()[w].name
                               << "  rule " << convert::to(share[w], sm3d)
                               << "   equations " << convert::to(r.well_rate[w], sm3d));
            rule_total += share[w];
            solved_total += r.well_rate[w];
        }
        BOOST_TEST_MESSAGE(what << ": target " << convert::to(target, sm3d)
                           << ", rule " << convert::to(rule_total, sm3d)
                           << ", equations " << convert::to(solved_total, sm3d));
        // The two must place the rate the same way. They need not reach the
        // target: a group asked for more than its wells can deliver gets what
        // they can, and both should say so rather than pretend.
        //
        // The cases with a well on a *rate* limit are reported rather than
        // asserted: they disagree, for the reason isolated in
        // a_rate_limited_well_stays_under_its_limit below.
        if (required) {
            for (std::size_t w = 0; w < share.size(); ++w) {
                BOOST_CHECK_CLOSE(convert::to(r.well_rate[w], sm3d),
                                  convert::to(share[w], sm3d), 1.0);
            }
            BOOST_CHECK_CLOSE(convert::to(solved_total, sm3d),
                              convert::to(rule_total, sm3d), 0.1);
        }
        BOOST_CHECK_LE(convert::to(solved_total, sm3d), convert::to(target, sm3d) * 1.001);
    };

    auto caseWithTarget = [&](const double fraction) {
        auto c = gnetinjeGas();
        double free_total = 0.0;
        for (const auto& w : c.wells()) {
            free_total += w.q_ref;
        }
        c.setGroupTarget(fraction * free_total);
        return std::make_pair(std::move(c), fraction * free_total);
    };

    // Every well able to take its share: the multiplier alone reproduces a
    // plain guide-rate split.
    {
        auto [c, target] = caseWithTarget(0.8);
        c.finish();
        compare("plain split", c, target, /*required=*/true);
    }

    // Hard against the target: a fifth of what the wells would take.
    {
        auto [c, target] = caseWithTarget(0.2);
        c.finish();
        compare("strongly binding", c, target, /*required=*/true);
    }

    // Barely binding -- the shares and the free rates almost coincide, which is
    // where any hysteresis in a control test will chatter.
    {
        auto [c, target] = caseWithTarget(0.999);
        c.finish();
        compare("marginally binding", c, target, /*required=*/true);
    }

    // One well that cannot take its share, so the rule has to redistribute.
    {
        auto [c, target] = caseWithTarget(1.0);
        c.wells()[0].rate_limit = convert::from(2.0e5, cubic(meter) / day);
        c.finish();
        compare("one well limited", c, target, /*required=*/false);
    }

    // Two of the four, on different branches, so the redistribution has to
    // cross the network as well as the group.
    {
        auto [c, target] = caseWithTarget(1.0);
        c.wells()[0].rate_limit = convert::from(2.0e5, cubic(meter) / day);
        c.wells()[2].rate_limit = convert::from(1.5e5, cubic(meter) / day);
        c.finish();
        compare("two wells limited", c, target, /*required=*/false);
    }

    // A target nobody can meet. The equations must not pretend otherwise.
    {
        auto [c, target] = caseWithTarget(2.0);
        c.finish();
        compare("beyond capacity", c, target, /*required=*/true);
    }
}

// How the formulations degrade as the wells stiffen. dq/dbhp sets the loop gain.
// Measured over the whole grid of starts, because a single start says too little.
BOOST_AUTO_TEST_CASE(stiffness_sweep)
{
    const auto n = static_cast<int>(startingPoints().size());
    for (const double stiffness : {1.0e4, 6.0e4, 3.0e5, 1.0e6}) {
        auto c = gnetinjeGas();
        c.setStiffness(stiffness);
        c.finish();
        BOOST_TEST_MESSAGE("dq/dbhp = " << stiffness << " sm3/d/bar");

        const EliminatedProblem eliminated_problem{c};
        const int bracket = basin("  bracketing (shipped)",
                                  [&](const State& p) { return bracketing(eliminated_problem, p, 0.1); });
        const int eliminated = basin("  eliminated, trust region",
                                     [&](const State& p) {
                                         return newton(eliminated_problem, p, TrustRegion{});
                                     });
        const int full = basin("  full, plain newton + bounds",
                               [&](const State& p) {
                                   FullProblem problem{c};
                                   problem.setEnforceBounds(true);
                                   return newton(problem, p, FullStep{});
                               });
        BOOST_CHECK_EQUAL(bracket, n);
        BOOST_CHECK_GE(eliminated, n - 1);
        BOOST_CHECK_GT(full, 3 * n / 4);
    }
}

// Trace one dumped system iteration by iteration: the rate every control offers
// each well, and the multiplier the group equation is currently implying. Set
// OPM_NETWORK_TRACE to one file written by --network-dump-failures. A cycling
// active set is unreadable from the outside and obvious from this.
BOOST_AUTO_TEST_CASE(trace_one_dumped_system)
{
    const char* path = std::getenv("OPM_NETWORK_TRACE");
    if (path == nullptr || !std::filesystem::is_regular_file(path)) {
        BOOST_TEST_MESSAGE("OPM_NETWORK_TRACE not set to a dump file, nothing to trace");
        return;
    }

    const auto gas = gnetinjeGas();
    std::ifstream in(path);
    auto [system, guess] = gas.systemFromDump(in);
    // solve() does this once before its loop; a trace that skips it is tracing
    // a different system.
    system.refreshGuides(system.start(guess));

    constexpr double perDay = 86400.0;
    constexpr double toBar = 1.0e-5;

    auto x = system.start(guess);
    for (int it = 1; it <= 16; ++it) {
        const bool moved = system.updateControls(x);
        const auto r = system.residual(x);
        double worst = 0.0;
        for (const auto e : r) {
            worst = std::max(worst, std::abs(e));
        }
        const double lambda = x[system.lambdaIdx()];

        std::ostringstream out;
        out << "it " << std::setw(2) << it << "  lambda " << std::setw(10) << lambda
            << "  |r| " << std::setw(10) << worst << (moved ? "   <- controls moved" : "");
        for (int w = 0; w < system.numWells(); ++w) {
            const auto& well = system.wells()[w];
            const double p = (well.node == 0) ? system.terminalPressure()
                                              : x[system.pIdx(well.node)];
            out << "\n      " << std::setw(5) << well.name
                << " [" << system.controlLetter(w) << "]"
                << "  p_node " << std::setw(7) << p * toBar
                << "  q "      << std::setw(9) << x[system.qwIdx(w)] * perDay
                << " | allows: thp " << std::setw(9) << system.thpPotential(well, p, /*cap_by_rate_limit=*/true) * perDay
                << "  bhp "          << std::setw(9) << (well.ipr_b * well.bhp_limit - well.ipr_a) * perDay
                << "  rate "         << std::setw(9) << well.rate_limit * perDay
                << "  grup "         << std::setw(9) << well.guide * lambda * perDay
                << "  (guide " << std::setw(9) << well.guide * perDay << ")";
        }
        BOOST_TEST_MESSAGE(out.str());

        if (worst < 1e-2 && !moved) {
            BOOST_TEST_MESSAGE("converged");
            return;
        }
        auto J = system.jacobian(x);
        std::vector<double> negative(system.size()), dx;
        for (int i = 0; i < system.size(); ++i) {
            negative[i] = -r[i];
        }
        BOOST_REQUIRE(J.solve(negative, dx));
        dx = system.limitStep(x, dx);
        for (int i = 0; i < system.size(); ++i) {
            x[i] += dx[i];
        }
    }
    BOOST_TEST_MESSAGE("did not settle in 16 iterations");
}

// A well on a rate limit stays under it, which took two goes.
//
// thpPotential() searches between the table's first rate and the smaller of the
// well's rate limit and the table's reach. That cap makes thp's allowance tie
// with the rate limit, and "smallest allowance wins" breaks the tie towards thp
// -- whose row is bhp = tableBhp(p_node, q) and says nothing about a rate, so
// the well settles wherever the tubing crosses the ipr, here 485 550 against a
// limit of 200 000.
//
// Simply removing the cap fixes this case and costs 511/529 of the globalisation
// basin down to 271/529: while the pressures are still moving a well's crossing
// routinely lies past its limit, rate control pins it there, and four wells
// pinned at their limits ask the network for several times what it carries.
// What works is to keep the cap while the solve is moving and drop it once there
// is a converged iterate to enforce the limit from -- see solve().
BOOST_AUTO_TEST_CASE(a_rate_limited_well_stays_under_its_limit)
{
    const auto sm3d = cubic(meter) / day;

    auto c = gnetinjeGas();
    double target = 0.0;
    for (const auto& w : c.wells()) {
        target += w.q_ref;
    }
    const double limit = convert::from(2.0e5, cubic(meter) / day);
    c.wells()[0].rate_limit = limit;
    c.setGroupTarget(target);
    c.finish();

    auto system = c.system();
    const auto r = NetworkSolve::solve(system, c.nodePressures(kStart), kParams, NetworkSolve::FullStep{});
    BOOST_REQUIRE(r.converged);
    BOOST_TEST_MESSAGE("rate-limited well: q " << convert::to(r.well_rate[0], sm3d)
                       << " against a limit of " << convert::to(limit, sm3d)
                       << ", control " << system.controlLetter(0));
    BOOST_CHECK_LE(convert::to(r.well_rate[0], sm3d), convert::to(limit, sm3d) * 1.001);
}

// Group control on the production network, the same two equations as injection:
// sum of the group's oil rates meets the target, and each held well takes
// guide * multiplier.
BOOST_AUTO_TEST_CASE(production_group_target_is_an_equation)
{
    const auto sm3d = cubic(meter) / day;

    ProductionCase c;
    const double free_total = c.freeTotal();
    const double target = 0.6 * free_total;
    c.setGroupTarget(target);

    auto system = c.system();
    const auto r = NetworkSolve::solve(system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
    BOOST_TEST_MESSAGE("production group: " << (r.converged ? "converged in " : "FAILED after ")
                       << r.iterations << " iterations, residual " << r.residual);
    BOOST_REQUIRE(r.converged);

    const double total = r.well_rate[0] + r.well_rate[1];
    BOOST_TEST_MESSAGE("free " << convert::to(free_total, sm3d) << ", target "
                       << convert::to(target, sm3d) << ", delivered " << convert::to(total, sm3d)
                       << " sm3/d  (" << convert::to(r.well_rate[0], sm3d) << " + "
                       << convert::to(r.well_rate[1], sm3d) << ")");
    BOOST_CHECK_CLOSE(convert::to(total, sm3d), convert::to(target, sm3d), 0.1);

    // Both wells held, so the split follows the guide rates.
    const auto& wells = system.wells();
    BOOST_CHECK_CLOSE(r.well_rate[0] / wells[0].guide, r.well_rate[1] / wells[1].guide, 0.1);
}

// The production group equations against the rule they replace, exactly the
// check the injection side gets: cap each well by what it can take at the
// pressures the equations settled on, share by guide rate, fix any well whose
// share exceeds its cap and re-divide.
BOOST_AUTO_TEST_CASE(production_equations_match_the_rule_based_allocation)
{
    const auto sm3d = cubic(meter) / day;

    auto ruleBased = [](const NetworkSolve::ProductionSystem<double>& system,
                        const std::vector<double>& node_pressure, const double target) {
        const auto& wells = system.wells();
        const int n = static_cast<int>(wells.size());
        std::vector<double> cap(n), share(n, 0.0);
        std::vector<bool> pooled(n, true);
        for (int w = 0; w < n; ++w) {
            const double p = node_pressure[wells[w].node];
            cap[w] = std::min({system.thpPotential(wells[w], p),
                               NetworkSolve::ProductionSystem<double>::ipr(wells[w], 1,
                                                                           wells[w].bhp_limit),
                               wells[w].oil_rate_limit});
        }
        double remaining = target;
        for (int pass = 0; pass <= n; ++pass) {
            double guides = 0.0;
            for (int w = 0; w < n; ++w) {
                if (pooled[w]) {
                    guides += wells[w].guide;
                }
            }
            if (guides <= 0.0) {
                break;
            }
            bool fixed_one = false;
            for (int w = 0; w < n; ++w) {
                if (!pooled[w]) {
                    continue;
                }
                const double sh = remaining * wells[w].guide / guides;
                if (sh > cap[w]) {
                    share[w] = cap[w];
                    pooled[w] = false;
                    remaining -= cap[w];
                    fixed_one = true;
                    break;
                }
                share[w] = sh;
            }
            if (!fixed_one) {
                break;
            }
        }
        return share;
    };

    ProductionCase base;
    const double free_total = base.freeTotal();

    for (const double fraction : {0.9, 0.6, 0.25, 0.999, 1.5}) {
        ProductionCase c;
        c.setGroupTarget(fraction * free_total);
        auto system = c.system();
        const auto r = NetworkSolve::solve(system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
        if (!r.converged) {
            BOOST_TEST_MESSAGE("fraction " << fraction << ": does not converge");
            BOOST_CHECK(false);
            continue;
        }
        const auto share = ruleBased(system, r.node_pressure, c.groupTarget());
        double rule_total = 0.0, solved_total = 0.0;
        for (std::size_t w = 0; w < share.size(); ++w) {
            rule_total += share[w];
            solved_total += r.well_rate[w];
        }
        BOOST_TEST_MESSAGE("fraction " << fraction << ": target "
                           << convert::to(c.groupTarget(), sm3d) << ", rule "
                           << convert::to(rule_total, sm3d) << ", equations "
                           << convert::to(solved_total, sm3d) << "  ("
                           << convert::to(r.well_rate[0], sm3d) << " + "
                           << convert::to(r.well_rate[1], sm3d) << ")");
        for (std::size_t w = 0; w < share.size(); ++w) {
            BOOST_CHECK_CLOSE(convert::to(r.well_rate[w], sm3d), convert::to(share[w], sm3d), 1.0);
        }
        BOOST_CHECK_LE(convert::to(solved_total, sm3d),
                       convert::to(c.groupTarget(), sm3d) * 1.001);
    }
}

// How far from the answer the production solve can start and still land on it,
// and what happens when a limit binds -- the production counterpart of
// globalisation_basin, and the measure any change to its control rule has to
// answer to. Four configurations over a wide grid of starting node pressures.
BOOST_AUTO_TEST_CASE(production_control_rule_basin)
{
    const auto sm3d = cubic(meter) / day;
    ProductionCase base;
    const double free_total = base.freeTotal();

    struct Config { const char* what; double bhp_limit; double oil_limit; double fraction; };
    const std::vector<Config> configs{
        {"free",         40.0, 600.0, 0.0},
        {"bhp binds",    85.0, 600.0, 0.0},
        {"rate binds",   40.0,  80.0, 0.0},
        {"group binds",  40.0, 600.0, 0.5},
    };

    int all_solved = 0, all_total = 0;
    for (const auto& cfg : configs) {
        int solved = 0, total = 0, iterations = 0;
        for (int pf = 0; pf < 14; ++pf) {
            const double p = convert::from(50.0 + 30.0 * pf, bars);
            ProductionCase c;
            for (auto& w : c.wells()) {
                w.bhp_limit = convert::from(cfg.bhp_limit, bars);
                w.oil_rate_limit = convert::from(cfg.oil_limit, sm3d);
            }
            c.setGroupTarget(cfg.fraction * free_total);
            auto system = c.system();
            const auto r = NetworkSolve::solve(
                system, std::vector<double>{convert::from(80.0, bars), p}, kParams, NetworkSolve::FullStep{});
            ++total;
            if (r.converged) {
                ++solved;
                iterations += r.iterations;
            }
        }
        BOOST_TEST_MESSAGE("production basin, " << std::setw(11) << cfg.what << "  "
                           << solved << "/" << total << " starts, mean "
                           << (solved > 0 ? iterations / solved : 0) << " iterations");
        all_solved += solved;
        all_total += total;
    }
    BOOST_TEST_MESSAGE("production basin, total       " << all_solved << "/" << all_total);
    BOOST_CHECK_GT(all_solved, 3 * all_total / 4);
}

// Does resolving the split actually flip less, or only converge more?
//
// The two selections on the same systems, over the whole 529-point grid of
// starting pressures at four group targets, counting the iterations on which
// some well changed control -- not just whether the solve landed. Taking the
// share from the iterate's multiplier is the rule that cycles; resolving it is
// the shipped one.
//
// The guides matter more than anything else here. Left at each well's own
// reference rate they are proportional to what each well can take, every share
// is feasible, nobody is ever dropped from the pool and the active set is
// uniformly GGGG -- there is nothing for either rule to get wrong, and the two
// measure the same. That is not the situation a timestep produces: guides are
// set once, from the previous operating point, and the solve then settles
// somewhere else. `equal_guides` is that case -- four equal guides against wells
// whose capacities differ by a factor of two, which is what the gas deck's own
// failures look like.
BOOST_AUTO_TEST_CASE(resolving_the_split_flips_less)
{
    const auto starts = startingPoints();

    struct Tally { int solved = 0; int switches = 0; int cycled = 0; };
    auto run = [&](const bool from_multiplier, const double fraction, const bool equal_guides) {
        Tally t;
        for (const auto& start : starts) {
            auto c = gnetinjeGas();
            double free_total = 0.0;
            for (const auto& w : c.wells()) {
                free_total += w.q_ref;
            }
            if (equal_guides) {
                for (auto& w : c.wells()) {
                    w.guide = free_total / c.wells().size();
                }
            }
            c.setGroupTarget(fraction * free_total);
            c.finish();
            auto system = c.system();
            system.setGroupShareFromMultiplier(from_multiplier);
            // nodePressures() spreads the two applied pressures over the whole
            // tree; System::start() wants one per node, and handing it the bare
            // pair reads past the end of it.
            const auto r = NetworkSolve::solve(system, c.nodePressures(start), kParams, NetworkSolve::FullStep{});
            t.switches += r.switches;
            if (r.converged) {
                ++t.solved;
            } else if (r.controls_moving) {
                ++t.cycled;
            }
        }
        return t;
    };

    for (const bool equal_guides : {false, true}) {
        const std::string heading = equal_guides ? "guides equal, capacities differ:"
                                                 : "guides proportional to capacity:";
        BOOST_TEST_MESSAGE(heading);
        int old_solved = 0, new_solved = 0;
        int old_switches = 0, new_switches = 0;
        int old_cycled = 0, new_cycled = 0;
        for (const double fraction : {0.7, 0.95, 1.0, 1.05, 1.2}) {
            const auto from_lambda = run(true, fraction, equal_guides);
            const auto resolved = run(false, fraction, equal_guides);
            BOOST_TEST_MESSAGE(
                "  target " << fraction << " of free:  multiplier "
                << from_lambda.solved << "/" << starts.size() << ", "
                << from_lambda.switches << " switches, " << from_lambda.cycled << " cycling"
                << "   |   resolved " << resolved.solved << "/" << starts.size() << ", "
                << resolved.switches << " switches, " << resolved.cycled << " cycling");
            old_solved += from_lambda.solved;   new_solved += resolved.solved;
            old_switches += from_lambda.switches; new_switches += resolved.switches;
            old_cycled += from_lambda.cycled;   new_cycled += resolved.cycled;
        }
        BOOST_TEST_MESSAGE("  totals: multiplier " << old_solved << " solved, " << old_switches
                           << " switches, " << old_cycled << " cycling  |  resolved "
                           << new_solved << " solved, " << new_switches << " switches, "
                           << new_cycled << " cycling");
        if (equal_guides) {
            // Where the guides do not already encode each well's capacity --
            // which is every real timestep -- all three have to improve.
            BOOST_CHECK_GT(new_solved, old_solved);
            BOOST_CHECK_LT(new_switches, old_switches);
            BOOST_CHECK_LT(new_cycled, old_cycled);
        }
    }
}

// One network the simulator could not solve, kept verbatim.
//
// Written by --network-dump-failures from GNETINJE_GAS-01 at a step where the
// field target sits just above what the four injectors can deliver. It is the
// shape every one of those failures had: equal guide rates against wells whose
// capacities differ by two, and a target the pool cannot meet. Taking the group
// share from the iterate's multiplier cycles TTGG / TTTT here and never lands;
// resolving the split converges in a handful of iterations.
BOOST_AUTO_TEST_CASE(a_dumped_simulator_failure_converges)
{
    const std::string dump = R"(phase GAS
terminal 3.4e+07
group_target 18.3248
guides_from_potential 1
analytic_jacobian 1
node PLAT-A -1 9999
node M5S 0 3
node G1 1 9999
node M5N 1 2
node F1 3 9999
well F-1H 4 1 39469 0.0013367 4.25e+07 11.5741 11.5741 4.96921 1
well F-2H 4 1 85174 0.00288533 4.25e+07 11.5741 11.5741 4.97161 1
well G-3H 2 1 69029.9 0.00233706 4.25e+07 11.5741 11.5741 4.19187 1
well G-4H 2 1 76010.4 0.00257402 4.25e+07 11.5741 11.5741 4.19211 1
guess 3.4e+07 4.50559e+07 4.50559e+07 4.994e+07 4.994e+07
)";

    const auto gas = gnetinjeGas();
    auto solve = [&](const bool from_multiplier) {
        std::istringstream in(dump);
        auto [system, guess] = gas.systemFromDump(in);
        system.setGroupShareFromMultiplier(from_multiplier);
        return NetworkSolve::solve(system, guess, kParams, NetworkSolve::FullStep{});
    };

    const auto from_lambda = solve(true);
    BOOST_TEST_MESSAGE("share from the multiplier: "
                       << (from_lambda.converged ? "converged in " : "FAILED after ")
                       << from_lambda.iterations << " iterations, " << from_lambda.switches
                       << " switches, controls " << from_lambda.control_trace);
    BOOST_CHECK(!from_lambda.converged);
    // It is still moving the set on two thirds of its iterations when it gives
    // up; `controls_moving` only reports the last one, which lands either way.
    BOOST_CHECK_GT(from_lambda.switches, from_lambda.iterations / 2);

    const auto resolved = solve(false);
    BOOST_TEST_MESSAGE("split resolved:            "
                       << (resolved.converged ? "converged in " : "FAILED after ")
                       << resolved.iterations << " iterations, " << resolved.switches
                       << " switches");
    BOOST_CHECK(resolved.converged);
    BOOST_CHECK_LT(resolved.iterations, 15);
    BOOST_CHECK_LT(resolved.switches, 5);
}

// Is the relaxed update ever needed as a fallback?
//
// Over the grid of starting pressures, ungrouped and at three group targets,
// with the step limiter the simulator runs -- and with a line search behind the
// full step, to see whether a globalisation would do instead.
//
// Everything up to a target the wells can just about meet is solved outright,
// and the line search never gets a chance to help. What is left, at 1.2 times
// what the wells can deliver, is not a globalisation problem: the retry recovers
// none of it. Those are an active set chasing the pressures -- everyone on group
// control drives the nodes out to where no well has a potential, every well is
// dropped from the pool, nobody is on group control, the pressures recover and
// the shares look feasible again: GGGG, TTTT, GGTT, round again. A different
// cycle from the one resolving the split removed, and with no multiplier in it.
//
// So the answer is: not on either deck -- both run with no fallback at all --
// but yes in general, and this is the case it is still there for.
BOOST_AUTO_TEST_CASE(the_fallback_still_has_one_case_to_cover)
{
    const auto starts = startingPoints();

    auto measure = [&](const double fraction) {
        int plain = 0, retried = 0, cycling = 0;
        for (const auto& start : starts) {
            auto c = gnetinjeGas();
            double free_total = 0.0;
            for (const auto& w : c.wells()) {
                free_total += w.q_ref;
            }
            c.setGroupTarget(fraction * free_total);
            c.finish();

            auto system = c.system();
            const auto guess = c.nodePressures(start);
            const auto first = NetworkSolve::solve(system, guess, kParams, NetworkSolve::FullStep{});
            if (first.converged) {
                ++plain;
                ++retried;
                continue;
            }
            auto again = c.system();
            const auto second =
                NetworkSolve::solve(again, guess, {1e-2, 50}, NetworkSolve::LineSearch{});
            if (second.converged) {
                ++retried;
            } else if (second.switches > second.iterations / 4) {
                ++cycling;
            }
        }
        BOOST_TEST_MESSAGE("target " << fraction << " of free:  full step " << plain
                           << "/" << starts.size() << ",  line search behind it "
                           << retried << "/" << starts.size() << ",  of the rest "
                           << cycling << " cycling");
        return std::make_tuple(plain, retried, cycling);
    };

    const auto n = static_cast<int>(starts.size());
    for (const double fraction : {0.0, 0.95, 1.05}) {
        const auto [plain, retried, cycling] = measure(fraction);
        BOOST_CHECK_EQUAL(plain, n);
        BOOST_CHECK_EQUAL(retried, n);
        BOOST_CHECK_EQUAL(cycling, 0);
    }

    // Over capacity: a line search buys nothing, and every failure is cycling
    // rather than a step that overshot.
    const auto [plain, retried, cycling] = measure(1.2);
    BOOST_CHECK_EQUAL(plain, retried);
    BOOST_CHECK_LT(retried, n);
    BOOST_CHECK_EQUAL(cycling, n - retried);
}

// thp that does not hold a well back must not be the control it ends on.
//
// With the bhp limit above what the tubing needs at the node pressure (and
// below shut-in, so the well still produces), the bhp limit binds.
// thpPotential() used to report exactly what the bhp limit allows, which
// tied, and the tie went to thp -- whose row then settled the bhp below the
// limit the deck set.
BOOST_AUTO_TEST_CASE(production_thp_that_does_not_bind_leaves_the_well_on_bhp)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    ProductionCase c;
    for (auto& w : c.wells()) {
        w.bhp_limit = convert::from(110.0, bars);
    }
    auto system = c.system();
    const auto r = NetworkSolve::solve(system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
    BOOST_REQUIRE(r.converged);
    for (int w = 0; w < system.numWells(); ++w) {
        const auto& well = system.wells()[w];
        BOOST_TEST_MESSAGE(well.name << " [" << system.controlLetter(w) << "] oil "
                           << convert::to(r.well_rate[w], cubic(meter) / day)
                           << ", bhp limit allows "
                           << convert::to(Sys::ipr(well, 1, well.bhp_limit), cubic(meter) / day));
        BOOST_CHECK_EQUAL(system.controlLetter(w), 'B');
        BOOST_CHECK_CLOSE(r.well_rate[w], Sys::ipr(well, 1, well.bhp_limit), 1e-6);
    }
}

// Efficiency factors, lift gas and node sources all change what the branch
// carries without changing what any well does. The check is the same for each:
// the node pressure is the table's answer to the branch flow the terms imply,
// and the wells' own rates are what they were without the term.
BOOST_AUTO_TEST_CASE(what_enters_a_branch_besides_the_wells)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    const auto sm3d = cubic(meter) / day;

    // Reference: nothing but the wells.
    ProductionCase plain;
    auto ref_system = plain.system();
    const auto ref = NetworkSolve::solve(ref_system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
    BOOST_REQUIRE(ref.converged);

    auto nodePressureFromBranch = [&](const Sys& system, const std::array<double, 3>& q) {
        return system.tableBhp(3, convert::from(80.0, bars), q, 0.0);
    };

    // Branch flow the wells deliver at a solution, from their bhp.
    auto wellTriple = [&](const Sys& system, const NetworkSolve::Result<double>& r, const int w) {
        // back out bhp from the oil rate, then the other phases from it
        const auto& well = system.wells()[w];
        const double bhp = (r.well_rate[w] - well.ipr_a[1]) / well.ipr_b[1];
        return std::array<double, 3>{Sys::ipr(well, 0, bhp), Sys::ipr(well, 1, bhp), Sys::ipr(well, 2, bhp)};
    };

    // 1. Well efficiency 0.5 on both wells: the branch carries half, the node
    //    pressure falls, the wells' own rates are whatever the new node
    //    pressure gives them -- and the node pressure is consistent with the
    //    halved branch.
    {
        ProductionCase c;
        for (auto& w : c.wells()) {
            w.efficiency = 0.5;
        }
        auto system = c.system();
        const auto r = NetworkSolve::solve(system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
        BOOST_REQUIRE(r.converged);
        std::array<double, 3> branch{};
        for (int w = 0; w < system.numWells(); ++w) {
            const auto q = wellTriple(system, r, w);
            for (int ph = 0; ph < 3; ++ph) {
                branch[ph] += 0.5 * q[ph];
            }
        }
        BOOST_TEST_MESSAGE("efficiency 0.5: node " << convert::to(r.node_pressure[1], bars)
                           << " bar against " << convert::to(ref.node_pressure[1], bars)
                           << " with the wells at full weight");
        BOOST_CHECK_LT(r.node_pressure[1], ref.node_pressure[1]);
        BOOST_CHECK_CLOSE(r.node_pressure[1], nodePressureFromBranch(system, branch), 1e-2);
    }

    // 2. Lift gas on one well: only the gas stream in the branch grows.
    {
        ProductionCase c;
        const double lift = convert::from(20000.0, sm3d);
        c.wells()[0].lift_gas = lift;
        auto system = c.system();
        const auto r = NetworkSolve::solve(system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
        BOOST_REQUIRE(r.converged);
        std::array<double, 3> branch{};
        for (int w = 0; w < system.numWells(); ++w) {
            const auto q = wellTriple(system, r, w);
            for (int ph = 0; ph < 3; ++ph) {
                branch[ph] += q[ph];
            }
        }
        branch[2] += lift;
        BOOST_TEST_MESSAGE("lift gas: node " << convert::to(r.node_pressure[1], bars) << " bar");
        BOOST_CHECK_CLOSE(r.node_pressure[1], nodePressureFromBranch(system, branch), 1e-2);
        // This table has one GFR point, so more gas cannot move its pressure.
        // Check the stream itself: the branch's gas at the starting point is
        // the wells' plus the lift, and start() is built the way residual() is.
        const auto with = system.start(ProductionCase::guess());
        const auto without = ref_system.start(ProductionCase::guess());
        BOOST_CHECK_CLOSE(with[system.qIdx(1, 2)] - without[ref_system.qIdx(1, 2)], lift, 1e-9);
    }

    // 3. A satellite at the node: a constant triple with no well behind it.
    {
        ProductionCase c;
        const std::array<double, 3> sat{convert::from(30.0, sm3d), convert::from(100.0, sm3d),
                                        convert::from(8000.0, sm3d)};
        auto system = c.system();
        system.setNodeSource(1, sat);
        const auto r = NetworkSolve::solve(system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
        BOOST_REQUIRE(r.converged);
        std::array<double, 3> branch = sat;
        for (int w = 0; w < system.numWells(); ++w) {
            const auto q = wellTriple(system, r, w);
            for (int ph = 0; ph < 3; ++ph) {
                branch[ph] += q[ph];
            }
        }
        BOOST_TEST_MESSAGE("satellite: node " << convert::to(r.node_pressure[1], bars) << " bar");
        BOOST_CHECK_CLOSE(r.node_pressure[1], nodePressureFromBranch(system, branch), 1e-2);
        BOOST_CHECK_GT(r.node_pressure[1], ref.node_pressure[1]);
    }
}

// The same for an injection network: a well at half efficiency halves what the
// branch carries, and the analytic Jacobian has to know it.
BOOST_AUTO_TEST_CASE(injection_efficiency_enters_the_branch)
{
    auto c = gnetinjeGas();
    for (auto& w : c.wells()) {
        w.efficiency = 0.5;
    }
    c.finish();
    auto system = c.system();
    system.setAnalyticJacobian(true);

    // A missed factor shows as a disagreement between the assembled Jacobian
    // and the differenced one, on the balance rows.
    const auto x = system.start(c.nodePressures(kStart));
    const auto J = system.jacobian(x);
    const auto r0 = system.residual(x);
    double worst = 0.0;
    for (int j = 0; j < system.size(); ++j) {
        auto shifted = x;
        const double h = 1e-3 * system.columnScale(j);
        shifted[j] += h;
        const auto rj = system.residual(shifted);
        for (int i = 0; i < system.size(); ++i) {
            worst = std::max(worst, std::abs(J(i, j) - (rj[i] - r0[i]) / h));
        }
    }
    BOOST_TEST_MESSAGE("analytic vs differenced with efficiencies: largest difference " << worst);
    BOOST_CHECK_LT(worst, 1e-3);

    // And it still solves, to a lower node pressure than at full weight.
    const auto r = NetworkSolve::solve(system, c.nodePressures(kStart), kParams, NetworkSolve::FullStep{});
    BOOST_REQUIRE(r.converged);
    auto full = gnetinjeGas();
    full.finish();
    auto full_system = full.system();
    const auto rf = NetworkSolve::solve(full_system, full.nodePressures(kStart), kParams, NetworkSolve::FullStep{});
    BOOST_REQUIRE(rf.converged);
    BOOST_TEST_MESSAGE("M5S at half efficiency " << convert::to(r.node_pressure[1], bars)
                       << " bar, at full " << convert::to(rf.node_pressure[1], bars));
    BOOST_CHECK_GT(r.node_pressure[1], rf.node_pressure[1]);
}

// An autochoke is a valve upstream of a node that throttles the oil collected
// there to a target; the node pressure is whatever that takes. Closed when the
// wells could deliver more than the target at the manifold pressure, open when
// they cannot -- and then the node is at the manifold pressure and the group
// produces what it can.
BOOST_AUTO_TEST_CASE(an_autochoke_holds_the_target_or_opens)
{
    const auto sm3d = cubic(meter) / day;
    ProductionCase base;
    const double free_total = base.freeTotal();

    // Closed: target below what the wells deliver freely.
    {
        ProductionCase c;
        auto system = c.system();
        system.setChokeTarget(1, 0.5 * free_total);
        const auto r = NetworkSolve::solve(system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
        BOOST_REQUIRE(r.converged);
        const double total = r.well_rate[0] + r.well_rate[1];
        BOOST_TEST_MESSAGE("choked: node " << convert::to(r.node_pressure[1], bars)
                           << " bar, oil " << convert::to(total, sm3d) << " against target "
                           << convert::to(0.5 * free_total, sm3d));
        BOOST_CHECK(system.choked(1));
        BOOST_CHECK_CLOSE(total, 0.5 * free_total, 0.1);
        BOOST_CHECK_GT(r.node_pressure[1], convert::from(80.0, bars));
        for (int w = 0; w < system.numWells(); ++w) {
            BOOST_CHECK_EQUAL(system.controlLetter(w), 'T');
        }
    }

    // Open: target beyond reach.
    {
        ProductionCase c;
        auto system = c.system();
        system.setChokeTarget(1, 2.0 * free_total);
        const auto r = NetworkSolve::solve(system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
        BOOST_REQUIRE(r.converged);
        const double total = r.well_rate[0] + r.well_rate[1];
        BOOST_TEST_MESSAGE("open: node " << convert::to(r.node_pressure[1], bars)
                           << " bar, oil " << convert::to(total, sm3d));
        BOOST_CHECK(!system.choked(1));
        BOOST_CHECK_CLOSE(total, free_total, 0.1);
    }
}

// The production Jacobian, assembled, against the differenced one -- on a plain
// network, under a group target, and with a closed choke, since each adds rows
// of its own. The rate derivatives come from differentiating the table lookup
// itself; a missed chain-rule term or a sign shows here and nowhere else.
BOOST_AUTO_TEST_CASE(production_analytic_jacobian_matches_differences)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    ProductionCase base;
    const double free_total = base.freeTotal();

    auto check = [&](const char* what, Sys& system) {
        system.setAnalyticJacobian(true);
        const auto x0 = system.start(ProductionCase::guess());
        // at the solution, so the active set the Jacobian is built on is the real one
        const auto r = NetworkSolve::solve(system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
        BOOST_REQUIRE(r.converged);
        auto x = x0;
        for (int n = 1; n <= system.numNodes(); ++n) { x[system.pIdx(n)] = r.node_pressure[n]; }
        system.updateControls(x);
        const auto J = system.jacobian(x);
        const auto r0 = system.residual(x);
        double worst = 0.0, largest = 0.0;
        for (int j = 0; j < system.size(); ++j) {
            auto shifted = x;
            const double h = 1e-4 * system.columnScale(j);
            shifted[j] += h;
            const auto rj = system.residual(shifted);
            for (int i = 0; i < system.size(); ++i) {
                const double fd = (rj[i] - r0[i]) / h;
                worst = std::max(worst, std::abs(J(i, j) - fd));
                largest = std::max(largest, std::abs(fd));
            }
        }
        BOOST_TEST_MESSAGE(what << ": largest entry " << largest << ", largest difference " << worst);
        BOOST_CHECK_LT(worst, 1e-3 * std::max(largest, 1.0));
    };

    { ProductionCase c; auto s = c.system(); check("plain", s); }
    { ProductionCase c; c.setGroupTarget(0.6 * free_total); auto s = c.system(); check("group target", s); }
    { ProductionCase c; auto s = c.system(); s.setChokeTarget(1, 0.5 * free_total); check("closed choke", s); }
}

// A well placed exactly where its tubing passes its own rate limit satisfies
// both rows, and an active set that picks one by the sign at the iterate flips
// between them. The tie goes to a Fischer-Burmeister row instead; it has to
// converge, land on the limit, and leave the tubing slack non-negative.
BOOST_AUTO_TEST_CASE(a_well_on_a_tie_converges)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    const auto sm3d = cubic(meter) / day;
    ProductionCase free_case;
    auto free_system = free_case.system();
    const auto free = NetworkSolve::solve(free_system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
    BOOST_REQUIRE(free.converged);

    // the limit is the rate the well freely produces: the tie, to the digit
    ProductionCase c;
    c.wells()[0].oil_rate_limit = free.well_rate[0];
    auto system = c.system();
    system.setAnalyticJacobian(true);
    const auto r = NetworkSolve::solve(system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
    BOOST_TEST_MESSAGE("tied well: " << (r.converged ? "converged in " : "FAILED after ")
                       << r.iterations << " iterations, " << r.switches << " switches, control '"
                       << system.controlLetter(0) << "', oil " << convert::to(r.well_rate[0], sm3d)
                       << " against limit " << convert::to(free.well_rate[0], sm3d));
    BOOST_CHECK(r.converged);
    BOOST_CHECK_LE(r.well_rate[0], free.well_rate[0] * (1.0 + 1e-3));
    BOOST_CHECK_CLOSE(r.well_rate[0], free.well_rate[0], 0.5);
}

// Do the assembled and the differenced Jacobian reach the same solution?
//
// Not "do both converge" -- the same node pressures and well rates, from the
// same start, over a grid of starting pressures, on every production system
// shape there is: plain, under a group target, with a closed choke, and with
// a well placed on its rate-limit tie. Where they part, one of them has found
// another root, and the simulator cannot tell which.
BOOST_AUTO_TEST_CASE(production_jacobians_reach_the_same_solution)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    ProductionCase base;
    const double free_total = base.freeTotal();
    auto free_system = base.system();
    const auto free = NetworkSolve::solve(free_system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
    BOOST_REQUIRE(free.converged);

    // The cases own the tables the systems point at, so they live here, not
    // inside the lambdas that build systems from them.
    ProductionCase plain_case, group_case, choke_case, tie_case;
    group_case.setGroupTarget(0.6 * free_total);
    tie_case.wells()[0].oil_rate_limit = free.well_rate[0];
    struct Shape { const char* what; std::function<Sys(void)> make; };
    const std::vector<Shape> shapes{
        {"plain",        [&] { return plain_case.system(); }},
        {"group target", [&] { return group_case.system(); }},
        {"closed choke", [&] { auto s = choke_case.system(); s.setChokeTarget(1, 0.5 * free_total); return s; }},
        {"on a tie",     [&] { return tie_case.system(); }},
    };

    for (const auto& shape : shapes) {
        int both = 0, agree = 0, only_fd = 0, only_an = 0, neither = 0;
        double worst_p = 0.0, worst_q = 0.0;
        for (int pf = 0; pf < 14; ++pf) {
            const std::vector<double> guess{convert::from(80.0, bars), convert::from(50.0 + 30.0 * pf, bars)};
            auto fd = shape.make(); fd.setAnalyticJacobian(false);
            auto an = shape.make(); an.setAnalyticJacobian(true);
            const auto rf = NetworkSolve::solve(fd, guess, kParams, NetworkSolve::FullStep{});
            const auto ra = NetworkSolve::solve(an, guess, kParams, NetworkSolve::FullStep{});
            if (rf.converged && ra.converged) {
                ++both;
                double dp = 0.0, dq = 0.0;
                for (std::size_t n = 0; n < rf.node_pressure.size(); ++n) {
                    dp = std::max(dp, std::abs(rf.node_pressure[n] - ra.node_pressure[n]));
                }
                for (std::size_t w = 0; w < rf.well_rate.size(); ++w) {
                    dq = std::max(dq, std::abs(rf.well_rate[w] - ra.well_rate[w])
                                      / std::max(std::abs(rf.well_rate[w]), 1e-12));
                }
                worst_p = std::max(worst_p, dp);
                worst_q = std::max(worst_q, dq);
                if (dp < convert::from(0.05, bars) && dq < 1e-3) { ++agree; }
            } else if (rf.converged) { ++only_fd; }
            else if (ra.converged) { ++only_an; }
            else { ++neither; }
        }
        BOOST_TEST_MESSAGE(std::setw(13) << shape.what << ": both " << both << "/14, agree " << agree
                           << ", only differenced " << only_fd << ", only analytic " << only_an
                           << ", neither " << neither << "; worst gap " << convert::to(worst_p, bars)
                           << " bar, " << 100 * worst_q << " % rate");
        BOOST_CHECK_EQUAL(agree, both);
        BOOST_CHECK_EQUAL(only_fd + only_an, 0);
    }
}

// Replay production systems the simulator wrote out (--network-dump-failures),
// against the MODEL5 tables. OPM_NETWORK_DUMP_PROD names the directory,
// OPM_VFP_INCLUDE the directory holding well_vfp.ecl and flowl_{b,c}_vfp.ecl.
// Each system is solved with both Jacobians; this is where a failure seen in
// the simulator becomes a millisecond of work.
BOOST_AUTO_TEST_CASE(replay_production_failures)
{
    const char* dir = std::getenv("OPM_NETWORK_DUMP_PROD");
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    if (dir == nullptr || inc == nullptr || !std::filesystem::is_directory(dir)) {
        BOOST_TEST_MESSAGE("OPM_NETWORK_DUMP_PROD / OPM_VFP_INCLUDE not set, nothing to replay");
        return;
    }
    std::deque<VFPProdTable> tables;
    VFPProdProperties<double> props;
    const UnitSystem units{};
    for (const char* name : {"well_vfp.ecl", "flowl_b_vfp.ecl", "flowl_c_vfp.ecl"}) {
        const auto path = std::filesystem::path(inc) / name;
        if (!std::filesystem::exists(path)) { continue; }
        const auto deck = Parser{}.parseFile(path.string());
        for (const auto& kw : deck.getKeywordList("VFPPROD")) {
            tables.emplace_back(*kw, /*gaslift_opt_active=*/true, units);
            props.addTable(tables.back());
        }
    }
    BOOST_TEST_MESSAGE("tables loaded: " << tables.size());

    std::vector<std::filesystem::path> files;
    for (const auto& e : std::filesystem::directory_iterator(dir)) {
        if (e.path().extension() == ".txt") { files.push_back(e.path()); }
    }
    std::sort(files.begin(), files.end());
    int fd_ok = 0, an_ok = 0, agree = 0, cm_ok = 0, cm_agree = 0;
    for (const auto& file : files) {
        std::ifstream in(file);
        std::string head; std::getline(in, head);
        if (head != "production") { continue; }
        auto [fd, guess] = NetworkSolve::readProduction<double>(in, props, units);
        auto an = fd; fd.setAnalyticJacobian(false); an.setAnalyticJacobian(true);
        auto cm = an; cm.setComplementarity(true);
        // OPM_NETWORK_MAX_IT raises the iteration cap, to tell "slow" from "stuck".
        const int max_it = std::getenv("OPM_NETWORK_MAX_IT") ? std::atoi(std::getenv("OPM_NETWORK_MAX_IT")) : 50;
        const auto rf = NetworkSolve::solve(fd, guess, {1e-2, max_it}, NetworkSolve::FullStep{});
        const auto ra = NetworkSolve::solve(an, guess, {1e-2, max_it}, NetworkSolve::FullStep{});
        const bool cm_ls = std::getenv("OPM_NETWORK_CM_LINESEARCH") != nullptr;
        const auto rc = cm_ls ? NetworkSolve::solve(cm, guess, {1e-2, max_it}, NetworkSolve::LineSearch{})
                              : NetworkSolve::solve(cm, guess, {1e-2, max_it}, NetworkSolve::FullStep{});
        fd_ok += rf.converged; an_ok += ra.converged; cm_ok += rc.converged;
        double gap = 0.0, cgap = 0.0;
        if (rf.converged && ra.converged) {
            for (std::size_t n = 0; n < rf.node_pressure.size(); ++n) {
                gap = std::max(gap, std::abs(rf.node_pressure[n] - ra.node_pressure[n]));
            }
            agree += (gap < convert::from(0.05, bars));
        }
        std::string where;
        if (ra.converged && rc.converged) {
            for (std::size_t n = 0; n < ra.node_pressure.size(); ++n) {
                cgap = std::max(cgap, std::abs(ra.node_pressure[n] - rc.node_pressure[n]));
            }
            cm_agree += (cgap < convert::from(0.05, bars));
            if (cgap >= convert::from(0.05, bars)) {
                for (int n = 1; n <= an.numNodes(); ++n) {
                    if (an.isChoke(n)) {
                        const int up = an.nodes()[n].parent;
                        where += fmt::format(" choke {} an:{} p {:.2f} (up {:.2f}) cm:{} p {:.2f}", an.nodes()[n].name,
                                             an.choked(n) ? "CLOSED" : "OPEN", ra.node_pressure[n] * 1e-5,
                                             (up == 0 ? an.terminalPressure() : ra.node_pressure[up]) * 1e-5,
                                             cm.choked(n) ? "CLOSED" : "OPEN", rc.node_pressure[n] * 1e-5);
                    }
                }
                for (int w = 0; w < an.numWells(); ++w) {
                    where += fmt::format(" {}:{}{:.0f}/{}{:.0f}", an.wells()[w].name, an.controlLetter(w),
                                         ra.well_rate[w] * 86400.0, cm.controlLetter(w), rc.well_rate[w] * 86400.0);
                }
            }
        }
        BOOST_TEST_MESSAGE("  " << file.filename().string()
                           << ": differenced " << (rf.converged ? "ok" : "FAILED") << " (" << rf.iterations << " it)"
                           << ", analytic " << (ra.converged ? "ok" : "FAILED") << " (" << ra.iterations << " it)"
                           << ", complementarity " << (rc.converged ? "ok" : "FAILED") << " (" << rc.iterations << " it)"
                           << (ra.converged && rc.converged ? fmt::format(", gap an/cm {:.3g} bar", convert::to(cgap, bars)) + where : "")
                           << (rc.converged ? "" : fmt::format("  cm residual {:.3g}", rc.residual)
                                               + (rc.control_trace.empty() ? "" : "  cm trace " + rc.control_trace)));
    }
    BOOST_TEST_MESSAGE("replayed " << files.size() << ": differenced ok " << fd_ok << ", analytic ok "
                       << an_ok << ", both ok and agreeing " << agree << "; complementarity ok " << cm_ok
                       << ", agreeing with analytic " << cm_agree);
}

// One complementarity row per well instead of an active set: on every shape,
// it must converge, and where the active set converges too the two must agree.
BOOST_AUTO_TEST_CASE(complementarity_agrees_with_the_active_set)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    const auto sm3d = cubic(meter) / day;
    ProductionCase base;
    const double free_total = base.freeTotal();
    auto free_system = base.system();
    const auto free = NetworkSolve::solve(free_system, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
    BOOST_REQUIRE(free.converged);

    ProductionCase plain_case, limited_case, choke_case, open_case, tie_case, bhp_case, high_case;
    for (auto& w : limited_case.wells()) { w.oil_rate_limit = 0.4 * free.well_rate[0]; }
    tie_case.wells()[0].oil_rate_limit = free.well_rate[0];
    for (auto& w : bhp_case.wells()) { w.bhp_limit = convert::from(110.0, bars); }
    for (auto& w : high_case.wells()) { w.bhp_limit = convert::from(150.0, bars); }   // above shut-in: nothing
    struct Shape { const char* what; std::function<Sys(void)> make; };
    const std::vector<Shape> shapes{
        {"plain",        [&] { return plain_case.system(); }},
        {"rate-limited", [&] { return limited_case.system(); }},
        {"bhp-limited",  [&] { return bhp_case.system(); }},
        {"bhp above reservoir", [&] { return high_case.system(); }},
        {"closed choke", [&] { auto s = choke_case.system(); s.setChokeTarget(1, 0.5 * free_total); return s; }},
        {"open choke",   [&] { auto s = open_case.system(); s.setChokeTarget(1, 2.0 * free_total); return s; }},
        {"on a tie",     [&] { return tie_case.system(); }},
    };
    for (const auto& shape : shapes) {
        int both = 0, agree = 0, only_as = 0, only_cm = 0, cm_its = 0;
        double worst = 0.0;
        for (int pf = 0; pf < 14; ++pf) {
            const std::vector<double> guess{convert::from(80.0, bars), convert::from(50.0 + 30.0 * pf, bars)};
            auto as = shape.make(); as.setAnalyticJacobian(true);
            auto cm = shape.make(); cm.setAnalyticJacobian(true); cm.setComplementarity(true);
            const auto ra = NetworkSolve::solve(as, guess, kParams, NetworkSolve::FullStep{});
            const auto rc = NetworkSolve::solve(cm, guess, kParams, NetworkSolve::FullStep{});
            if (rc.converged) { cm_its += rc.iterations; }
            const bool choke_shape = std::string(shape.what) == "closed choke" || std::string(shape.what) == "open choke";
            const bool disagree = ra.converged && rc.converged
                && std::abs(ra.well_rate[0] - rc.well_rate[0]) > 1e-3 * std::abs(ra.well_rate[0]);
            if (pf % 4 == 0 && (!rc.converged || disagree)) {
                std::string rates;
                for (std::size_t w = 0; w < rc.well_rate.size(); ++w) {
                    rates += fmt::format(" {:.0f}({}, bhp {:.1f}/as {:.1f})", rc.well_rate[w] * 86400.0,
                                         static_cast<int>(cm.control(static_cast<int>(w))), rc.well_bhp[w] * 1e-5, ra.well_bhp[w] * 1e-5);
                }
                BOOST_TEST_MESSAGE("    " << shape.what << ", start " << 50 + 30 * pf << " bar: cm " << (rc.converged ? "ok after " : "FAILED after ")
                                   << rc.iterations << " it, residual " << rc.residual << ", choke "
                                   << (cm.choked(1) ? "CLOSED" : "OPEN") << " p " << rc.node_pressure[1] * 1e-5
                                   << ", rates" << rates << ", trace " << rc.control_trace
                                   << " | active set: p " << ra.node_pressure[1] * 1e-5 << " rates "
                                   << ra.well_rate[0] * 86400.0 << " " << ra.well_rate[1] * 86400.0);
            }
            if (ra.converged && rc.converged) {
                ++both;
                double d = 0.0;
                for (std::size_t w = 0; w < ra.well_rate.size(); ++w) {
                    d = std::max(d, std::abs(ra.well_rate[w] - rc.well_rate[w]) / std::max(std::abs(ra.well_rate[w]), 1e-12));
                }
                worst = std::max(worst, d);
                agree += (d < 1e-3);
            } else if (ra.converged) { ++only_as; }
            else if (rc.converged) { ++only_cm; }
        }
        BOOST_TEST_MESSAGE(std::setw(13) << shape.what << ": both " << both << "/14, agree " << agree
                           << ", only active-set " << only_as << ", only complementarity " << only_cm
                           << ", worst rate gap " << 100 * worst << " %, complementarity mean its "
                           << (both + only_cm > 0 ? cm_its / (both + only_cm) : 0));
        BOOST_CHECK_EQUAL(only_as, 0);
        BOOST_CHECK_EQUAL(agree, both);
    }
}

// The production step bounds, checked directly: a 1000 bar node step is cut to
// 50, and a step that would take a node below an atmosphere stops at it.
// A well the well model has at zero rate above some node pressure -- the
// adapter's dead_above -- must not be produced from. The active set hands it
// to bhp control and reports 2500 m3/d through tubing that cannot lift
// (dump_prod_1112: three of them under a 6000 choke target); the
// complementarity shuts it. Both shapes from the dumps, on the bench tables.
BOOST_AUTO_TEST_CASE(complementarity_shuts_dead_wells)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    ProductionCase plain_case;
    auto plain = plain_case.system(); plain.setAnalyticJacobian(true);
    const auto free = NetworkSolve::solve(plain, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
    BOOST_REQUIRE(free.converged);
    const double p_free = free.node_pressure[1];

    // One dying well: dead just below where the node settles with it alive.
    {
        ProductionCase c;
        c.wells()[0].dead_above = p_free - convert::from(0.5, bars);
        auto cm = c.system(); cm.setAnalyticJacobian(true); cm.setComplementarity(true);
        auto as = c.system(); as.setAnalyticJacobian(true);
        const auto rc = NetworkSolve::solve(cm, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
        const auto ra = NetworkSolve::solve(as, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
        BOOST_REQUIRE(rc.converged);
        BOOST_TEST_MESSAGE("dying well: cm " << cm.controlLetter(0) << " " << rc.well_rate[0] * 86400
                           << " m3/d at node " << rc.node_pressure[1] * 1e-5 << " bar; active set "
                           << (ra.converged ? as.controlLetter(0) : '?') << " "
                           << (ra.converged ? ra.well_rate[0] * 86400 : 0.0) << " m3/d");
        BOOST_CHECK_EQUAL(cm.controlLetter(0), 'S');
        BOOST_CHECK_SMALL(rc.well_rate[0], 1e-12);
        BOOST_CHECK_GT(rc.well_rate[1], 0.0);
        // The other well alone leaves the node below where well 0 died --
        // the no-fixed-point shape; the shut has to stick.
        BOOST_CHECK_LT(rc.node_pressure[1], c.wells()[0].dead_above);
    }
    // A choke over wells that are all dead: nothing to throttle, valve open,
    // no flow. The row must not chase a target nothing can meet.
    {
        ProductionCase c;
        for (auto& w : c.wells()) { w.dead_above = convert::from(1.0, bars); }
        auto cm = c.system(); cm.setAnalyticJacobian(true); cm.setComplementarity(true);
        cm.setChokeTarget(1, 0.5 * free.well_rate[0]);
        const auto rc = NetworkSolve::solve(cm, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
        BOOST_REQUIRE(rc.converged);
        for (int w = 0; w < cm.numWells(); ++w) {
            BOOST_CHECK_EQUAL(cm.controlLetter(w), 'S');
            BOOST_CHECK_SMALL(rc.well_rate[w], 1e-12);
        }
        BOOST_CHECK_LT(rc.iterations, 10);
        BOOST_TEST_MESSAGE("choke over dead wells: " << rc.iterations << " it, node "
                           << rc.node_pressure[1] * 1e-5 << " bar");
    }
}

// The default path's hole, pinned so a fix flips this test rather than
// passing unnoticed: the active set gives a dead well (dead_above below the
// node pressure, thp unavailable) bhp control and produces from it.
BOOST_AUTO_TEST_CASE(the_active_set_still_produces_from_a_dead_well)
{
    ProductionCase c;
    c.wells()[0].dead_above = convert::from(1.0, bars);
    auto as = c.system(); as.setAnalyticJacobian(true);
    const auto r = NetworkSolve::solve(as, ProductionCase::guess(), kParams, NetworkSolve::FullStep{});
    BOOST_REQUIRE(r.converged);
    BOOST_TEST_MESSAGE("dead well on the active set: " << as.controlLetter(0) << " "
                       << r.well_rate[0] * 86400 << " m3/d");
    BOOST_CHECK_EQUAL(as.controlLetter(0), 'B');     // should one day be 'S'
    BOOST_CHECK_GT(r.well_rate[0], 0.0);
}

// The four AUTOCHK dumps behind the 2026-08-23 fixes, kept under
// tests/network_dumps. They need the model5 tables (OPM_VFP_INCLUDE); without
// them the case reports and passes, like replay_production_failures.
BOOST_AUTO_TEST_CASE(the_dumps_behind_the_complementarity_fixes_converge)
{
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    if (inc == nullptr) {
        BOOST_TEST_MESSAGE("OPM_VFP_INCLUDE not set; dumps not replayed");
        return;
    }
    std::deque<VFPProdTable> tables;
    VFPProdProperties<double> props;
    const UnitSystem units{};
    for (const char* name : {"well_vfp.ecl", "flowl_b_vfp.ecl", "flowl_c_vfp.ecl"}) {
        const auto deck = Parser{}.parseFile((std::filesystem::path(inc) / name).string());
        for (const auto& kw : deck.getKeywordList("VFPPROD")) {
            tables.emplace_back(*kw, true, units);
            props.addTable(tables.back());
        }
    }
    const auto dir = std::filesystem::path(__FILE__).parent_path() / "network_dumps";
    for (const char* name : {"autochk_prod_992.txt", "autochk_prod_1112.txt", "autochk_prod_0.txt", "autochk_prod_302.txt"}) {
        std::ifstream in(dir / name);
        BOOST_REQUIRE(in);
        std::string head; std::getline(in, head);
        auto [system, guess] = NetworkSolve::readProduction<double>(in, props, units);
        system.setAnalyticJacobian(true); system.setComplementarity(true);
        const auto r = NetworkSolve::solve(system, guess, {1e-2, 50}, NetworkSolve::FullStep{});
        std::string controls;
        for (int w = 0; w < system.numWells(); ++w) { controls += system.controlLetter(w); }
        BOOST_TEST_MESSAGE(name << ": " << (r.converged ? "ok" : "FAILED") << " in " << r.iterations
                           << " it, controls " << controls);
        BOOST_CHECK(r.converged);
        BOOST_CHECK_LT(r.iterations, 12);
    }
}

// thpPotential() searches for the IPR/tubing crossing -- 96 samples in bhp
// then 40 bisections. It does not have to: with the phase fractions held
// constant the IPR is linear in FLO and the table piecewise-linear on its own
// flow axis, so the crossing is exact per interval, which is what
// VFPHelpers::intersectWithIPR does. This sweeps node pressure over both and
// reports where they differ and what each costs. Needs OPM_VFP_INCLUDE for the
// model5 tables; reports and returns without them.
BOOST_AUTO_TEST_CASE(exact_intersection_against_the_scan)
{
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    if (inc == nullptr) {
        BOOST_TEST_MESSAGE("OPM_VFP_INCLUDE not set; scan/exact comparison skipped");
        return;
    }
    std::deque<VFPProdTable> tables;
    VFPProdProperties<double> props;
    const UnitSystem units{};
    for (const char* name : {"well_vfp.ecl", "flowl_b_vfp.ecl", "flowl_c_vfp.ecl"}) {
        const auto deck = Parser{}.parseFile((std::filesystem::path(inc) / name).string());
        for (const auto& kw : deck.getKeywordList("VFPPROD")) {
            tables.emplace_back(*kw, true, units);
            props.addTable(tables.back());
        }
    }
    // A real well, from the AUTOCHK dumps: table 1, and an inflow that puts the
    // crossing inside the table rather than off its end.
    const auto dir = std::filesystem::path(__FILE__).parent_path() / "network_dumps";
    std::ifstream in(dir / "autochk_prod_992.txt");
    BOOST_REQUIRE(in);
    std::string head; std::getline(in, head);
    auto [system, guess] = NetworkSolve::readProduction<double>(in, props, units);

    int compared = 0, differ = 0;
    double worst = 0.0, worst_p = 0.0;
    long scan_lookups = 0, exact_lookups = 0;
    for (int w = 0; w < system.numWells(); ++w) {
        const auto& well = system.wells()[w];
        for (int i = 0; i <= 40; ++i) {
            const double p = convert::from(15.0 + 1.5 * i, bars);
            system.resetLookups();
            const double a = system.thpPotentialFor(well, p);
            scan_lookups += system.lookups();
            system.resetLookups();
            const double b = system.thpPotentialExactFor(well, p, well.bhp_limit);
            exact_lookups += system.lookups();
            const bool inf_a = a == std::numeric_limits<double>::max();
            const bool inf_b = b == std::numeric_limits<double>::max();
            ++compared;
            if (inf_a || inf_b) {
                if (inf_a != inf_b) { ++differ; }
                continue;
            }
            const double scale = std::max({std::abs(a), std::abs(b), 1e-9});
            const double rel = std::abs(a - b) / scale;
            // Where they disagree, refine the scan: if it walks toward the
            // exact value the scan's resolution was the error, and if it does
            // not the difference is the constant-fraction assumption.
            if (!inf_a && !inf_b && std::abs(a - b) > 1e-3 * std::max({std::abs(a), std::abs(b), 1e-9})) {
                BOOST_TEST_MESSAGE("    " << well.name << " at " << convert::to(p, bars)
                                   << " bar: scan " << a * 86400.0 << ", exact " << b * 86400.0
                                   << " m3/d");
            }
            if (rel > 1e-3) {
                ++differ;
                if (rel > worst) { worst = rel; worst_p = convert::to(p, bars); }
            }
        }
    }
    BOOST_TEST_MESSAGE("compared " << compared << " (well, node pressure) pairs: "
                       << differ << " differ by more than 0.1 %"
                       << ", worst " << 100.0 * worst << " % at " << worst_p << " bar");
    // The scan's cost is table evaluations; the exact one's is axis walks, each
    // of which costs (flo axis + 1) evaluations inside VFPHelpers.
    const int axis = static_cast<int>(props.getTable(1).getFloAxis().size());
    const double exact_evals = double(exact_lookups) * (axis + 1);
    BOOST_TEST_MESSAGE("cost: scan " << scan_lookups << " table evaluations, exact "
                       << exact_lookups << " axis walks x (" << axis << "+1) = " << exact_evals
                       << "  (" << double(scan_lookups) / exact_evals << "x)");
    BOOST_CHECK_LT(exact_lookups, scan_lookups);
}

// How stale is the IPR by the time the network has finished with it?
//
// The IPR is an exact derivative of the converged well equation at one bhp --
// the point q_start sits at. Everything the network does with a well is a
// first-order prediction from there, so the honest question is how far it
// moves wells away from that point. If the moves are small the linearisation
// is fine and a richer inflow model would buy nothing; if they are large it is
// the limit on the whole approach.
//
// Reads a directory of dumps (OPM_NETWORK_DUMP_PROD) written with
// OPM_NETWORK_DUMP_ALL, so these are ordinary solves, not just failures.
BOOST_AUTO_TEST_CASE(how_far_the_network_moves_wells_from_their_ipr)
{
    const char* dir = std::getenv("OPM_NETWORK_DUMP_PROD");
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    if (dir == nullptr || inc == nullptr || !std::filesystem::is_directory(dir)) {
        BOOST_TEST_MESSAGE("OPM_NETWORK_DUMP_PROD / OPM_VFP_INCLUDE not set; staleness not measured");
        return;
    }
    std::deque<VFPProdTable> tables;
    VFPProdProperties<double> props;
    const UnitSystem units{};
    for (const char* name : {"well_vfp.ecl", "flowl_b_vfp.ecl", "flowl_c_vfp.ecl"}) {
        const auto path = std::filesystem::path(inc) / name;
        if (!std::filesystem::exists(path)) { continue; }
        const auto deck = Parser{}.parseFile(path.string());
        for (const auto& kw : deck.getKeywordList("VFPPROD")) {
            tables.emplace_back(*kw, true, units);
            props.addTable(tables.back());
        }
    }
    std::vector<std::filesystem::path> files;
    for (const auto& e : std::filesystem::directory_iterator(dir)) {
        if (e.path().extension() == ".txt") { files.push_back(e.path()); }
    }
    std::sort(files.begin(), files.end());

    std::vector<double> rate_move, bhp_move_bar, span_move;
    int systems = 0, wells = 0, from_zero = 0;
    for (const auto& file : files) {
        std::ifstream f(file);
        std::string head; std::getline(f, head);
        if (head != "production") { continue; }
        auto [system, guess] = NetworkSolve::readProduction<double>(f, props, units);
        system.setAnalyticJacobian(true);
        system.setComplementarity(true);
        const auto r = NetworkSolve::solve(system, guess, {1e-2, 50}, NetworkSolve::FullStep{});
        if (!r.converged) { continue; }
        ++systems;
        for (int w = 0; w < system.numWells(); ++w) {
            const auto& well = system.wells()[w];
            if (!(well.ipr_b[1] < 0.0)) { continue; }
            ++wells;
            if (!(well.q_start > 0.0)) { ++from_zero; continue; }
            const double q0 = well.q_start;
            const double q1 = r.well_rate[w];
            rate_move.push_back(std::abs(q1 - q0) / q0);
            const double bhp0 = (q0 - well.ipr_a[1]) / well.ipr_b[1];
            const double dbhp = std::abs(r.well_bhp[w] - bhp0);
            bhp_move_bar.push_back(dbhp * 1e-5);
            // Against the well's whole operating range, bhp limit to shut-in.
            // (Dividing by the drawdown from q_start instead just reproduces
            // |dq|/q_start exactly, the IPR being linear -- no new information.)
            double shut = well.bhp_limit;
            for (int ph = 0; ph < 3; ++ph) {
                if (well.ipr_b[ph] < 0.0) {
                    shut = std::max(shut, -well.ipr_a[ph] / well.ipr_b[ph]);
                }
            }
            if (shut > well.bhp_limit) {
                span_move.push_back(dbhp / (shut - well.bhp_limit));
            }
        }
    }
    auto pct = [](std::vector<double> v, const double q) {
        if (v.empty()) { return 0.0; }
        std::sort(v.begin(), v.end());
        return v[std::min(v.size() - 1, std::size_t(q * v.size()))];
    };
    BOOST_TEST_MESSAGE("solved " << systems << " systems, " << wells << " wells ("
                       << from_zero << " starting from zero rate, skipped)");
    BOOST_TEST_MESSAGE("  |dq|/q_start   median " << 100 * pct(rate_move, 0.5)
                       << " %, p90 " << 100 * pct(rate_move, 0.9)
                       << " %, p99 " << 100 * pct(rate_move, 0.99)
                       << " %, max " << 100 * pct(rate_move, 1.0) << " %");
    BOOST_TEST_MESSAGE("  |dbhp|/range   median " << 100 * pct(span_move, 0.5)
                       << " %, p90 " << 100 * pct(span_move, 0.9)
                       << " %, p99 " << 100 * pct(span_move, 0.99)
                       << " %, max " << 100 * pct(span_move, 1.0) << " %");
    BOOST_TEST_MESSAGE("  |dbhp|         median " << pct(bhp_move_bar, 0.5)
                       << " bar, p90 " << pct(bhp_move_bar, 0.9)
                       << " bar, p99 " << pct(bhp_move_bar, 0.99)
                       << " bar, max " << pct(bhp_move_bar, 1.0) << " bar");
}

BOOST_AUTO_TEST_CASE(the_group_tree_as_sparse_rows)
{
    TreeCase tc;
    const std::vector<double> guess{convert::from(120.0, bars)};

    struct Want { const char* what; double field; double p1; std::array<double,3> caps;
                  std::array<double,3> q;
                  TreeCase::Sys::Mode mode = TreeCase::Sys::Mode::Oil;
                  double wct = 0.0; };
    const std::vector<Want> cases{
        // Field 3000 shared 1:1 by the platforms; P2 can only make 400, so P1
        // absorbs 2600 and splits it evenly between two 2000-capable wells.
        {"one well caps, field target met", 3000.0, 0.0, {2000, 2000, 400}, {1300, 1300, 400}},
        // Now P1 is itself capped at 2000. Field cannot reach 3000: 2000 + 400.
        {"platform target binds too",       3000.0, 2000.0, {2000, 2000, 400}, {1000, 1000, 400}},
        // Nothing binds below the field: 1200 split 600/600, then 300/300.
        {"nothing below the field binds",   1200.0, 0.0, {2000, 2000, 2000}, {300, 300, 600}},
        // Same physical answer as the first case, asked for as liquid: every
        // well makes 25 % of its oil as water, so LRAT 3750 == ORAT 3000.
        {"a liquid target, 25 % water",     3750.0, 0.0, {2000, 2000, 400}, {1300, 1300, 400},
         TreeCase::Sys::Mode::Liquid, 0.25},
    };

    for (const auto& c : cases) {
        auto sys = tc.build(c.field, c.p1, c.caps, c.mode, c.wct);
        const auto r = NetworkSolve::solve(sys, guess, {1e-2, 60}, NetworkSolve::FullStep{});
        std::string got;
        for (int w = 0; w < sys.numWells(); ++w) {
            got += fmt::format(" {}={:.1f}", sys.wells()[w].name, r.well_rate[w] * 86400.0);
        }
        BOOST_TEST_MESSAGE(c.what << ": " << (r.converged ? "converged" : "FAILED")
                           << " in " << r.iterations << " it, residual " << r.residual << got);
        BOOST_CHECK(r.converged);
        if (!r.converged) { continue; }
        for (int w = 0; w < sys.numWells(); ++w) {
            BOOST_CHECK_CLOSE(r.well_rate[w] * 86400.0, c.q[w], 0.5);
        }
    }
}

// The same four shapes, resolved by an active set instead of by
// Fischer-Burmeister.
//
// Both paths solve the same equations; they differ only in who decides which
// limit holds each group and each well. With the choice made outside the
// Newton every group row is linear -- the definition rows always were, and the
// allocation row becomes one of three equalities -- so the group block is
// exact and constant. The point of the comparison is that the two agree on the
// answer; the iteration counts say what the switching costs.
//
// This is also the path that matters near convergence, where the active set is
// already right and there is nothing to switch.
BOOST_AUTO_TEST_CASE(the_group_tree_by_active_set)
{
    TreeCase tc;
    const std::vector<double> guess{convert::from(120.0, bars)};

    struct Want { const char* what; double field; double p1; std::array<double,3> caps;
                  std::array<double,3> q;
                  TreeCase::Sys::Mode mode = TreeCase::Sys::Mode::Oil;
                  double wct = 0.0; };
    const std::vector<Want> cases{
        {"one well caps, field target met", 3000.0, 0.0, {2000, 2000, 400}, {1300, 1300, 400}},
        {"platform target binds too",       3000.0, 2000.0, {2000, 2000, 400}, {1000, 1000, 400}},
        {"nothing below the field binds",   1200.0, 0.0, {2000, 2000, 2000}, {300, 300, 600}},
        {"a liquid target, 25 % water",     3750.0, 0.0, {2000, 2000, 400}, {1300, 1300, 400},
         TreeCase::Sys::Mode::Liquid, 0.25},
    };

    for (const auto& c : cases) {
        auto fb = tc.build(c.field, c.p1, c.caps, c.mode, c.wct, false);
        auto as = tc.build(c.field, c.p1, c.caps, c.mode, c.wct, true);
        const auto rf = NetworkSolve::solve(fb, guess, {1e-2, 60}, NetworkSolve::FullStep{});
        const auto ra = NetworkSolve::solve(as, guess, {1e-2, 60}, NetworkSolve::FullStep{});
        std::string got;
        for (int w = 0; w < as.numWells(); ++w) {
            got += fmt::format(" {}={:.1f}", as.wells()[w].name, ra.well_rate[w] * 86400.0);
        }
        BOOST_TEST_MESSAGE(c.what << ": active set "
                           << (ra.converged ? "converged" : "FAILED")
                           << " in " << ra.iterations << " it (fb "
                           << (rf.converged ? "converged" : "FAILED") << " in "
                           << rf.iterations << "), residual " << ra.residual << got);
        BOOST_CHECK(ra.converged);
        if (!ra.converged) { continue; }
        for (int w = 0; w < as.numWells(); ++w) {
            BOOST_CHECK_CLOSE(ra.well_rate[w] * 86400.0, c.q[w], 0.5);
            if (rf.converged) {
                BOOST_CHECK_CLOSE(ra.well_rate[w] * 86400.0, rf.well_rate[w] * 86400.0, 0.5);
            }
        }
    }
}

// The group tree built from a deck rather than by hand.
//
// tests/GROUPTREE.DATA: FIELD -> PLAT -> {GP1 -> W1 W2, GP2 -> W3}, PLAT on
// GCONPROD ORAT 3000, W3 held to 400 by its own limit. Reads GRUPTREE,
// GCONPROD and WCONPROD out of the Schedule and builds the sparse system from
// them, so the deck is the single description of the case -- and the same
// Schedule is what Stein's balancer needs for its GuideRate, which is what
// makes it usable as an oracle on exactly this configuration.
BOOST_AUTO_TEST_CASE(the_group_tree_from_a_deck)
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

    using Sys = NetworkSolve::ProductionSystem<double>;
    VFPProdProperties<double> props;
    const UnitSystem units{};
    Sys sys(props, units);

    NetworkSolve::Node term; term.name = "TERM"; term.parent = -1;
    sys.addNode(term, 0.0);
    NetworkSolve::Node n; n.name = "N"; n.parent = 0;
    sys.addNode(n, 0.0);
    sys.setTerminalPressure(convert::from(120.0, bars));

    // Groups, parents first so a parent's index is known when its child is
    // added. FIELD has no production target here; PLAT carries it.
    std::map<std::string, int> gidx;
    std::function<void(const std::string&, int)> addTree =
        [&](const std::string& name, const int parent) {
            typename Sys::Group g;
            g.name = name;
            g.parent = parent;
            const auto& grp = schedule.getGroup(name, 0);
            g.efficiency = grp.getGroupEfficiencyFactor();
            if (grp.isProductionGroup()) {
                const auto ctl = grp.productionControls(summary_state);
                if (ctl.cmode == Opm::Group::ProductionCMode::ORAT) {
                    g.mode = Sys::Mode::Oil;   g.target = ctl.oil_target;
                } else if (ctl.cmode == Opm::Group::ProductionCMode::LRAT) {
                    g.mode = Sys::Mode::Liquid; g.target = ctl.liquid_target;
                } else if (ctl.cmode == Opm::Group::ProductionCMode::GRAT) {
                    g.mode = Sys::Mode::Gas;   g.target = ctl.gas_target;
                } else if (ctl.cmode == Opm::Group::ProductionCMode::WRAT) {
                    g.mode = Sys::Mode::Water; g.target = ctl.water_target;
                }
            }
            const int me = sys.addGroup(std::move(g));
            gidx[name] = me;
            for (const auto& child : grp.groups()) { addTree(child, me); }
        };
    addTree("FIELD", -1);
    BOOST_TEST_MESSAGE("groups from the deck: " << gidx.size());

    // Wells, with their own ORAT ceiling from WCONPROD as the individual limit.
    for (const auto& wname : schedule.wellNames(0)) {
        const auto& well = schedule.getWell(wname, 0);   // Opm::Well
        if (!well.isProducer()) { continue; }
        const auto ctl = well.productionControls(summary_state);
        typename Sys::Well w;
        w.name = wname;
        w.node = 1;
        w.vfp_table = 0;                       // no tubing on this deck
        // A linear inflow reaching the well's own ORAT limit at 100 bar and
        // shutting in at 250; the deck fixes the limit, this fixes the slope.
        const double cap = ctl.oil_rate;
        const double b = -cap / (250e5 - 100e5);
        w.ipr_a[1] = -b * 250e5;
        w.ipr_b[1] = b;
        w.bhp_limit = convert::from(100.0, bars);
        w.oil_rate_limit = cap;
        w.guide = 1.0;
        w.group = gidx.at(well.groupName());
        w.q_start = 0.5 * cap;
        sys.addWell(std::move(w));
    }
    sys.setGroupTree(true);
    sys.setAnalyticJacobian(true);
    sys.setComplementarity(true);
    sys.finish();
    sys.finishGroups();

    const auto r = NetworkSolve::solve(sys, {convert::from(120.0, bars)},
                                       {1e-2, 60}, NetworkSolve::FullStep{});
    std::string got;
    double total = 0.0;
    for (int w = 0; w < sys.numWells(); ++w) {
        got += fmt::format(" {}={:.1f}", sys.wells()[w].name, r.well_rate[w] * 86400.0);
        total += r.well_rate[w] * 86400.0;
    }
    BOOST_TEST_MESSAGE("deck-driven: " << (r.converged ? "converged" : "FAILED")
                       << " in " << r.iterations << " it," << got
                       << ", PLAT total " << total);
    BOOST_REQUIRE(r.converged);
    // PLAT's ORAT target is 3000 and the wells can make 4400, so it binds.
    BOOST_CHECK_CLOSE(total, 3000.0, 0.5);
    // W3 is capped at 400 by its own limit, so GP1 carries the remaining 2600.
    for (int w = 0; w < sys.numWells(); ++w) {
        const double q = r.well_rate[w] * 86400.0;
        const double want = sys.wells()[w].name == "W3" ? 400.0 : 1300.0;
        BOOST_CHECK_CLOSE(q, want, 0.5);
    }
}

BOOST_AUTO_TEST_SUITE_END()
