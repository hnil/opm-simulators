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
#include <opm/simulators/wells/network/NetworkSteinStart.hpp>

// The decks and the dumps: trees built from a deck, generated instances,
// and every route scored by the judge on the systems MODEL5 dumped.

BOOST_AUTO_TEST_SUITE(NetworkDumpsBench)

// The 200-well instance where the two routes stop at different answers,
// each self-consistent by its own walk: the groups of both, the totals under
// every target, and whether the outer loop started from the reduced answer
// stays there.
BOOST_AUTO_TEST_CASE(generated_large_disagreement)
{
    using Sys = DeckTrees::Sys;
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    std::string what;
    const auto text = generateDeck({200, 40, 30, 2u}, &what);
    DeckTrees dt(text, DeckTrees::FromText{});
    for (const double j : {1.0, 2.5}) {
        DeckTrees::Ipr ipr; ipr.j_scale = j;
        auto reduced = dt.build(0, ipr);
        auto outer = dt.build(0, ipr);
        auto again = dt.build(0, ipr);
        const auto rr = NetworkSolve::solveReduced(reduced, dt.guess(20.0), params, true);
        const auto ro = NetworkSolve::solveWithTree(outer, dt.guess(20.0), params, NetworkSolve::FullStep{});
        BOOST_REQUIRE(rr.converged && ro.result.converged);
        const auto ra = NetworkSolve::solveWithTree(again, rr.node_pressure, params, NetworkSolve::FullStep{});
        const auto stein_r = steinAllocation(reduced, reduced.reducedState(), *dt.schedule, 0);
        const auto stein_o = steinAllocation(outer, ro.result.state, *dt.schedule, 0);
        reduced.updateControls(reduced.reducedState());
        outer.updateControls(ro.result.state);
        int differ = 0, again_differ = 0;
        double total_r = 0.0, total_o = 0.0;
        for (int w = 0; w < reduced.numWells(); ++w) {
            total_r += rr.well_rate[w]; total_o += ro.result.well_rate[w];
            if (std::abs(rr.well_rate[w] - ro.result.well_rate[w]) > 0.005 * std::max(rr.well_rate[w], 1e-9)) { ++differ; }
            if (ra.result.converged && std::abs(rr.well_rate[w] - ra.result.well_rate[w]) > 0.005 * std::max(rr.well_rate[w], 1e-9)) { ++again_differ; }
        }
        BOOST_TEST_MESSAGE(fmt::format(
            "j {}: {}; reduced {} it, outer {} passes / {} it; {} of {} wells differ; totals {:.0f} / {:.0f};"
            " outer from the reduced answer: {} passes / {} it{}, {} wells differ from reduced",
            j, what, rr.iterations, ro.passes, ro.inner_iterations, differ, reduced.numWells(),
            total_r * 86400.0, total_o * 86400.0, ra.passes, ra.inner_iterations,
            ra.consistent ? "" : " INCONSISTENT", again_differ));
        const auto& groups = reduced.groups();
        auto under = [&](const Sys& sys, const std::vector<double>& q, const int g) {
            double t = 0.0;
            for (int w = 0; w < sys.numWells(); ++w) {
                for (int a = sys.wells()[w].group; a >= 0; a = groups[a].parent) {
                    if (a == g) { t += q[w]; break; }
                }
            }
            return t;
        };
        std::vector<double> qs_r(stein_r.size()), qs_o(stein_o.size());
        for (std::size_t w = 0; w < stein_r.size(); ++w) { qs_r[w] = stein_r[w] / 86400.0; }
        for (std::size_t w = 0; w < stein_o.size(); ++w) { qs_o[w] = stein_o[w] / 86400.0; }
        std::string gl;
        for (int g = 0; g < reduced.numGroups(); ++g) {
            if (!(groups[g].target > 0.0)) { continue; }
            auto letter = [](const Sys::GroupBind b) { return b == Sys::GroupBind::Free ? 'F' : b == Sys::GroupBind::Own ? 'O' : 'S'; };
            gl += fmt::format("\n    {:5} <- {:5} {} {:7.0f}: reduced {:7.0f} [{}] stein {:7.0f} | outer {:7.0f} [{}] stein {:7.0f}",
                              groups[g].name, groups[g].parent >= 0 ? groups[groups[g].parent].name : "-",
                              groups[g].mode == Sys::Mode::Liquid ? "LRAT" : "ORAT", groups[g].target * 86400.0,
                              under(reduced, rr.well_rate, g) * 86400.0 * (groups[g].mode == Sys::Mode::Liquid ? 1.25 : 1.0),
                              letter(reduced.groupBind(g)),
                              stein_r.empty() ? 0.0 : under(reduced, qs_r, g) * 86400.0 * (groups[g].mode == Sys::Mode::Liquid ? 1.25 : 1.0),
                              under(outer, ro.result.well_rate, g) * 86400.0 * (groups[g].mode == Sys::Mode::Liquid ? 1.25 : 1.0),
                              letter(outer.groupBind(g)),
                              stein_o.empty() ? 0.0 : under(outer, qs_o, g) * 86400.0 * (groups[g].mode == Sys::Mode::Liquid ? 1.25 : 1.0));
        }
        BOOST_TEST_MESSAGE("  targets (liquid totals for LRAT, oil otherwise):" << gl);
        // Where the two differ: the wells' capacities at each answer's
        // pressures, and how many wells the tubing cannot lift at all there.
        auto anatomy = [&](const Sys& sys, const std::vector<double>& p, const std::vector<double>& q) {
            int dead = 0, held = 0, thp = 0, own = 0;
            double cap_sum = 0.0, p_max = 0.0;
            for (int w = 0; w < sys.numWells(); ++w) {
                const auto c = sys.control(w);
                held += c == Sys::Control::Tree ? 1 : 0;
                thp += c == Sys::Control::Thp ? 1 : 0;
                own += (c != Sys::Control::Tree && c != Sys::Control::Thp) ? 1 : 0;
                const double allow = sys.ownAllowance(w);
                cap_sum += allow < 1e30 ? allow : 0.0;
                const auto& well = sys.wells()[w];
                if (well.vfp_table > 0 && !(sys.thpPotential(well, p[well.node]) > 0.0)) { ++dead; }
            }
            for (std::size_t n = 1; n < p.size(); ++n) { p_max = std::max(p_max, convert::to(p[n], bars)); }
            double total = 0.0;
            for (const double x : q) { total += x; }
            return fmt::format("total {:.0f}, capacity sum {:.0f}, {} held / {} thp / {} own, {} dead, node p max {:.1f} bar",
                               total * 86400.0, cap_sum * 86400.0, held, thp, own, dead, p_max);
        };
        BOOST_TEST_MESSAGE("  reduced: " << anatomy(reduced, rr.node_pressure, rr.well_rate));
        BOOST_TEST_MESSAGE("  outer:   " << anatomy(outer, ro.result.node_pressure, ro.result.well_rate));
        // The same wells differ: which groups they sit under, and their capacities at both.
        std::map<std::string, int> where;
        int cap_bigger_at_reduced = 0, cap_bigger_at_outer = 0;
        for (int w = 0; w < reduced.numWells(); ++w) {
            if (std::abs(rr.well_rate[w] - ro.result.well_rate[w]) <= 0.005 * std::max(rr.well_rate[w], 1e-9)) { continue; }
            ++where[groups[reduced.wells()[w].group].name];
            const double cr = reduced.ownAllowance(w), co = outer.ownAllowance(w);
            if (cr > co * 1.005) { ++cap_bigger_at_reduced; } else if (co > cr * 1.005) { ++cap_bigger_at_outer; }
        }
        std::string wl;
        for (const auto& [g, k] : where) { wl += fmt::format(" {}:{}", g, k); }
        BOOST_TEST_MESSAGE(fmt::format("  differing wells by group:{}; capacity larger at the reduced answer for {}, at the outer for {}",
                                       wl, cap_bigger_at_reduced, cap_bigger_at_outer));
    }
}

// Step 4: real inflow performance. The dumps in tests/network_dumps are
// MODEL5 production systems the simulator wrote out -- every well's IPR,
// guide, limits and tubing datum as the well model had them, the choke on
// B1 holding 6000. Here the choke is turned off and the deck's group tree
// put on instead, so B1's 6000 is held by allocation. Two mechanisms for
// one target: the choke raises B1's pressure until the wells self-limit on
// their tubing; the tree assigns the rates and lets B1's pressure be what
// the branch says. Both are answers; they need not be the same one.
BOOST_AUTO_TEST_CASE(model5_dumps_with_the_tree)
{
    using Sys = DeckTrees::Sys;
    const std::string model5 = kNetworkDecks + "NETWORK_MODEL5_STDW_AUTOCHK.DATA";
    const auto dumps = std::filesystem::path(__FILE__).parent_path() / "network_dumps";
    if (!std::filesystem::exists(model5) || !std::filesystem::is_directory(dumps)) {
        BOOST_TEST_MESSAGE("opm-tests or the dumps not present, skipping");
        return;
    }
    // The tables the dumps refer to, with lift gas active as the dumps assume.
    std::deque<VFPProdTable> tables;
    VFPProdProperties<double> props;
    const UnitSystem units{};
    for (const char* name : {"well_vfp.ecl", "flowl_b_vfp.ecl", "flowl_c_vfp.ecl"}) {
        const auto path = std::filesystem::path(kNetworkDecks) / "include" / name;
        if (!std::filesystem::exists(path)) { continue; }
        const auto deck = Parser{}.parseFile(path.string());
        for (const auto& kw : deck.getKeywordList("VFPPROD")) {
            tables.emplace_back(*kw, /*gaslift_opt_active=*/true, units);
            props.addTable(tables.back());
        }
    }
    DeckTrees dt(model5);
    const int step = 1;
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    std::vector<std::filesystem::path> files;
    for (const auto& e : std::filesystem::directory_iterator(dumps)) {
        if (e.path().extension() == ".txt") { files.push_back(e.path()); }
    }
    std::sort(files.begin(), files.end());
    for (const auto& file : files) {
        std::ifstream in(file);
        std::string head; std::getline(in, head);
        if (head != "production") { continue; }
        auto [dumped, guess] = NetworkSolve::readProduction<double>(in, props, units);
        dumped.setAnalyticJacobian(true);
        dumped.setComplementarity(true);
        int b1 = -1;
        for (int n = 0; n <= dumped.numNodes(); ++n) { if (dumped.nodes()[n].name == "B1") { b1 = n; } }
        // As dumped: the choke holds B1.
        auto choke = dumped;
        const auto rc = NetworkSolve::solve(choke, guess, params, NetworkSolve::FullStep{});
        // The tree instead.
        auto tree = dumped;
        if (b1 >= 0) { tree.setChokeTarget(b1, 0.0); }
        dt.attachTree(tree, step);
        auto outer = tree, reduced = tree, fb = tree;
        const auto ro = NetworkSolve::solveWithTree(outer, guess, params, NetworkSolve::FullStep{});
        const auto rr = NetworkSolve::solveReduced(reduced, guess, params, true);
        const auto rf = NetworkSolve::solve(fb, guess, params, NetworkSolve::FullStep{});
        std::string sets;
        for (const auto& s : ro.sets) { sets += " " + s; }
        // What the walk says at the outer answer, when it disagrees with the set solved.
        std::string at_answer;
        if (ro.result.converged && !ro.consistent) {
            outer.updateControls(ro.result.state);
            outer.resolveTree(ro.result.state);
            at_answer = " -> at the answer " + outer.treeSignature();
        }
        BOOST_TEST_MESSAGE(fmt::format(
            "{}: {} wells, tree{}; choke-only {} in {} it; tree: outer {} passes / {} it{} [{}]{}, reduced {} it,"
            " fb {}; off the tables: {}/{} ({}/{})",
            file.filename().string(), dumped.numWells(), dt.description,
            rc.converged ? "converged" : "FAILED", rc.iterations, ro.passes, ro.inner_iterations,
            ro.consistent ? "" : " NOT CONSISTENT", sets, at_answer, rr.iterations,
            rf.converged ? std::to_string(rf.iterations) + " it" : "FAILED", ro.result.off_axis, rr.off_axis,
            outer.offAxisNote(), reduced.offAxisNote()));
        auto line = [&](const char* tag, const Sys& sys, const std::vector<double>& q,
                        const std::vector<double>& p, const std::vector<double>& bhp) {
            std::string l = fmt::format("  {:9}", tag);
            double under_b1 = 0.0;
            for (int w = 0; w < sys.numWells(); ++w) {
                l += fmt::format(" {}={:.0f}({}{:.0f})", sys.wells()[w].name, q[w] * 86400.0,
                                 sys.controlLetter(w), bhp.empty() ? 0.0 : convert::to(bhp[w], bars));
                if (sys.wells()[w].node == b1) { under_b1 += q[w]; }
            }
            l += "  nodes";
            for (int n = 1; n <= sys.numNodes(); ++n) {
                l += fmt::format(" {}={:.1f}", sys.nodes()[n].name, convert::to(p[n], bars));
            }
            l += fmt::format("  at B1 {:.0f}", under_b1 * 86400.0);
            BOOST_TEST_MESSAGE(l);
            return under_b1 * 86400.0;
        };
        const double at_b1_choke = rc.converged ? line("choke", choke, rc.well_rate, rc.node_pressure, rc.well_bhp) : 0.0;
        const double at_b1_tree = ro.result.converged ? line("tree", outer, ro.result.well_rate, ro.result.node_pressure, ro.result.well_bhp) : 0.0;
        if (ro.result.converged) {
            // A thp well's rate in the outer answer against the crossing the
            // reduced form would put it on at the same node pressure.
            std::string gaps;
            for (int w = 0; w < outer.numWells(); ++w) {
                if (outer.control(w) != Sys::Control::Thp) { continue; }
                const auto& well = outer.wells()[w];
                const double cross = outer.thpPotential(well, ro.result.node_pressure[well.node]);
                gaps += fmt::format(" {} newton {:.0f} crossing {:.0f}", well.name,
                                    ro.result.well_rate[w] * 86400.0, cross * 86400.0);
            }
            if (!gaps.empty()) { BOOST_TEST_MESSAGE("  thp wells:" << gaps); }
        }
        BOOST_CHECK(rc.converged);
        BOOST_CHECK(ro.result.converged);
        BOOST_CHECK(ro.consistent);
        BOOST_CHECK(rr.converged);
        if (!ro.result.converged || !rr.converged) { continue; }
        // An answer off the tables is reported, not compared: nothing says
        // what the tables mean there.
        if (ro.result.off_axis == 0 && rr.off_axis == 0) {
            for (int w = 0; w < outer.numWells(); ++w) {
                BOOST_CHECK_CLOSE(ro.result.well_rate[w], rr.well_rate[w], 0.5);
                if (rf.converged) { BOOST_WARN_CLOSE(ro.result.well_rate[w], rf.well_rate[w], 0.5); }
            }
        }
        const auto stein = steinAllocation(reduced, reduced.reducedState(), *dt.schedule, step);
        BOOST_CHECK(!stein.empty());
        if (!stein.empty()) {
            std::string both;
            double worst = 0.0;
            for (int w = 0; w < reduced.numWells(); ++w) {
                const double mine = rr.well_rate[w] * 86400.0;
                both += fmt::format(" {}={:.0f}/{:.0f}", reduced.wells()[w].name, mine, stein[w]);
                worst = std::max(worst, std::abs(mine - stein[w]) / std::max(std::abs(stein[w]), 1.0));
            }
            BOOST_TEST_MESSAGE("  reduced/stein:" << both);
            BOOST_CHECK_LT(worst, 0.005);
        }
        // Both mechanisms hold B1 at 6000 when the wells can make it.
        if (rc.converged && at_b1_choke > 5990.0) { BOOST_CHECK_CLOSE(at_b1_tree, 6000.0, 0.5); }
    }
}

// The simulator's rules as a fourth route, on the dumps the simulator
// itself wrote when those rules failed, and on the two decks.
BOOST_AUTO_TEST_CASE(the_legacy_rules_on_the_dumps)
{
    using Sys = DeckTrees::Sys;
    const std::string model5 = kNetworkDecks + "NETWORK_MODEL5_STDW_AUTOCHK.DATA";
    const auto dumps = std::filesystem::path(__FILE__).parent_path() / "network_dumps";
    if (!std::filesystem::exists(model5) || !std::filesystem::is_directory(dumps)) {
        BOOST_TEST_MESSAGE("opm-tests or the dumps not present, skipping");
        return;
    }
    std::deque<VFPProdTable> tables;
    VFPProdProperties<double> props;
    const UnitSystem units{};
    for (const char* name : {"well_vfp.ecl", "flowl_b_vfp.ecl", "flowl_c_vfp.ecl"}) {
        const auto path = std::filesystem::path(kNetworkDecks) / "include" / name;
        if (!std::filesystem::exists(path)) { continue; }
        const auto deck = Parser{}.parseFile(path.string());
        for (const auto& kw : deck.getKeywordList("VFPPROD")) {
            tables.emplace_back(*kw, /*gaslift_opt_active=*/true, units);
            props.addTable(tables.back());
        }
    }
    DeckTrees dt(model5);
    const int step = 1;
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    auto report = [&](const std::string& what, Sys& legacy, const std::vector<double>& guess,
                      Sys* tree_sys, const NetworkSolve::TreeResult<double>* tree) {
        const auto rl = NetworkSolve::solveLegacy(legacy, guess);
        std::string l = fmt::format("{}: legacy {} in {} it, imbalance {:.3f} bar, {} lookups off the tables", what,
                                    rl.converged ? "converged" : "NOT converged", rl.iterations,
                                    convert::to(rl.imbalance, bars), rl.off_axis);
        if (rl.choke_thp) {
            l += fmt::format("; choke thp {:.2f} bar ({}), mismatch {:+.1f} %, {} evaluations",
                             convert::to(*rl.choke_thp, bars), rl.choke_bracketed ? "bracketed" : "kept",
                             100 * rl.choke_mismatch, rl.choke_evaluations);
        } else if (rl.choke_evaluations > 0) {
            l += fmt::format("; choke open (no bracket in {} evaluations)", rl.choke_evaluations);
        }
        BOOST_TEST_MESSAGE(l);
        std::string w = "  legacy   ";
        double total = 0.0;
        for (int i = 0; i < legacy.numWells(); ++i) {
            w += fmt::format(" {}={:.0f}({})", legacy.wells()[i].name, rl.well_rate[i] * 86400.0, rl.controls[i]);
            total += rl.well_rate[i];
        }
        w += "  nodes";
        for (int n = 1; n <= legacy.numNodes(); ++n) {
            w += fmt::format(" {}={:.1f}", legacy.nodes()[n].name, convert::to(rl.node_pressure[n], bars));
        }
        w += fmt::format("  total {:.0f}", total * 86400.0);
        BOOST_TEST_MESSAGE(w);
        if (tree && tree->result.converged) {
            std::string t = "  tree     ";
            double tt = 0.0;
            for (int i = 0; i < tree_sys->numWells(); ++i) {
                t += fmt::format(" {}={:.0f}({})", tree_sys->wells()[i].name, tree->result.well_rate[i] * 86400.0,
                                 tree_sys->controlLetter(i));
                tt += tree->result.well_rate[i];
            }
            t += "  nodes";
            for (int n = 1; n <= tree_sys->numNodes(); ++n) {
                t += fmt::format(" {}={:.1f}", tree_sys->nodes()[n].name, convert::to(tree->result.node_pressure[n], bars));
            }
            t += fmt::format("  total {:.0f}", tt * 86400.0);
            BOOST_TEST_MESSAGE(t);
        }
        return rl;
    };

    std::vector<std::filesystem::path> files;
    for (const auto& e : std::filesystem::directory_iterator(dumps)) {
        if (e.path().extension() == ".txt") { files.push_back(e.path()); }
    }
    std::sort(files.begin(), files.end());
    for (const auto& file : files) {
        std::ifstream in(file);
        std::string head; std::getline(in, head);
        if (head != "production") { continue; }
        auto [dumped, guess] = NetworkSolve::readProduction<double>(in, props, units);
        auto legacy = dumped;
        auto tree = dumped;
        tree.setAnalyticJacobian(true);
        tree.setComplementarity(true);
        for (int n = 0; n <= tree.numNodes(); ++n) { if (tree.nodes()[n].name == "B1") { tree.setChokeTarget(n, 0.0); } }
        dt.attachTree(tree, step);
        const auto ro = NetworkSolve::solveWithTree(tree, guess, params, NetworkSolve::FullStep{});
        report(file.filename().string(), legacy, guess, &tree, &ro);
    }

    // The decks: MODEL5 with its choke on, NETWORK-01 with nothing binding.
    {
        auto sys = dt.build(step, {});
        for (int n = 0; n <= sys.numNodes(); ++n) { if (sys.nodes()[n].name == "B1") { sys.setChokeTarget(n, 6000.0 / 86400.0); } }
        auto tree = dt.build(step, {});
        const auto ro = NetworkSolve::solveWithTree(tree, dt.guess(21.0), params, NetworkSolve::FullStep{});
        const auto rl = report("MODEL5 deck, choke on", sys, dt.guess(21.0), &tree, &ro);
        BOOST_CHECK(rl.converged);
    }
    const std::string n01 = kNetworkDecks + "NETWORK-01.DATA";
    if (std::filesystem::exists(n01)) {
        DeckTrees d1(n01);
        auto sys = d1.build(step, {});
        auto tree = d1.build(step, {});
        const auto ro = NetworkSolve::solveWithTree(tree, d1.guess(80.0), params, NetworkSolve::FullStep{});
        const auto rl = report("NETWORK-01 deck", sys, d1.guess(80.0), &tree, &ro);
        BOOST_CHECK(rl.converged);
        BOOST_REQUIRE(ro.result.converged);
        for (int w = 0; w < sys.numWells(); ++w) {
            BOOST_CHECK_CLOSE(rl.well_rate[w], ro.result.well_rate[w], 1.0);
        }
        for (int n = 1; n <= sys.numNodes(); ++n) {
            BOOST_CHECK_CLOSE(rl.node_pressure[n], ro.result.node_pressure[n], 0.5);
        }
    }
}

// Every network system a simulator run wrote out (OPM_NETWORK_DUMP_ALL), by
// three routes: the full system's active set, the reduced form, the legacy
// rules. Where they disagree on whether a well flows, the disagreement is
// printed with the well's crossing at each answer's pressure -- the
// dead-or-alive multiplicity, if that is what it is.
BOOST_AUTO_TEST_CASE(replay_dumps_by_three_routes)
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
    std::vector<std::filesystem::path> files;
    for (const auto& e : std::filesystem::directory_iterator(dir)) {
        if (e.path().extension() == ".txt") { files.push_back(e.path()); }
    }
    auto number = [](const std::filesystem::path& p) {
        const auto s = p.stem().string();
        return std::atoi(s.substr(s.rfind('_') + 1).c_str());
    };
    std::sort(files.begin(), files.end(), [&](const auto& a, const auto& b) { return number(a) < number(b); });
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    int n = 0, full_ok = 0, red_ok = 0, leg_ok = 0, hold_ok = 0, held_wells = 0, dead_disagree = 0, shown = 0, shown_failed = 0,
        failed_grouped = 0, failed_plain = 0, shown_failed_grouped = 0;
    std::map<std::string, int> dead_full, dead_red, dead_leg, dead_hold;
    for (const auto& file : files) {
        std::ifstream in(file);
        std::string head; std::getline(in, head);
        if (head != "production") { continue; }
        auto [dumped, guess] = NetworkSolve::readProduction<double>(in, props, units);
        dumped.setAnalyticJacobian(true);
        dumped.setComplementarity(true);
        auto full = dumped, red = dumped, leg = dumped, hold = dumped;
        const auto rf = NetworkSolve::solve(full, guess, params, NetworkSolve::FullStep{});
        const auto rr = NetworkSolve::solveReduced(red, guess, params, true);
        const auto rl = NetworkSolve::solveLegacy(leg, guess);
        const auto rh = NetworkSolve::solveReduced(hold, guess, params, true, NetworkSolve::CliffRule::Hold);
        ++n; full_ok += rf.converged; red_ok += rr.converged; leg_ok += rl.converged;
        hold_ok += rh.converged; held_wells += rh.held_at_cliff;
        for (int w = 0; w < dumped.numWells(); ++w) {
            if (rh.converged && !(rh.well_rate[w] > 1e-9) && !dumped.wells()[w].shut) { dead_hold[dumped.wells()[w].name] += 1; }
        }
        if (!rr.converged) {
            (dumped.groupTarget() > 0.0 ? failed_grouped : failed_plain) += 1;
            if ((dumped.groupTarget() > 0.0 ? shown_failed_grouped : shown_failed)++ < 5) {
                std::string dead;
                for (int w = 0; w < dumped.numWells(); ++w) { dead += red.controlLetter(w); }
                BOOST_TEST_MESSAGE(fmt::format("  {} reduced NOT converged{}: {} it, {} stalls, {} cliffs, residual {:.3g}, set {} [{}]",
                                               file.filename().string(), dumped.groupTarget() > 0.0 ? " (flat target)" : "",
                                               rr.iterations, rr.stalls, rr.cliffs, rr.residual, dead, rr.sets.substr(0, 90)));
            }
        }
        bool disagree = false;
        std::string detail;
        for (int w = 0; w < dumped.numWells(); ++w) {
            const auto& well = dumped.wells()[w];
            if (well.shut) { continue; }
            const bool df = rf.converged && !(rf.well_rate[w] > 1e-9);
            const bool dr = rr.converged && !(rr.well_rate[w] > 1e-9);
            const bool dl = rl.converged && !(rl.well_rate[w] > 1e-9);
            dead_full[well.name] += df; dead_red[well.name] += dr; dead_leg[well.name] += dl;
            if (rf.converged && rr.converged && rl.converged && (df != dr || df != dl)) {
                disagree = true;
                const double pf = rf.node_pressure[well.node], pr = rr.node_pressure[well.node], pl = rl.node_pressure[well.node];
                detail += fmt::format("\n    {}: full {:.0f} at {:.2f} bar [{}], reduced {:.0f} at {:.2f} [{}], legacy {:.0f} at {:.2f} [{}];"
                                      " crossing at those pressures {:.0f} / {:.0f} / {:.0f}; dead_above {:.2f}, q_start {:.0f}",
                                      well.name, rf.well_rate[w] * 86400.0, convert::to(pf, bars), full.controlLetter(w),
                                      rr.well_rate[w] * 86400.0, convert::to(pr, bars), red.controlLetter(w),
                                      rl.well_rate[w] * 86400.0, convert::to(pl, bars), rl.controls[w],
                                      dumped.thpPotential(well, pf) * 86400.0, dumped.thpPotential(well, pr) * 86400.0,
                                      dumped.thpPotential(well, pl) * 86400.0,
                                      convert::to(well.dead_above, bars), well.q_start * 86400.0);
            }
        }
        if (disagree) {
            ++dead_disagree;
            if (shown++ < 6) { BOOST_TEST_MESSAGE(file.filename().string() << ":" << detail); }
        }
    }
    BOOST_TEST_MESSAGE(fmt::format("{} dumps: converged full {} / reduced {} / legacy {} / reduced-hold {} ({} wells held at a cliff);"
                                   " {} with a dead-or-alive disagreement; reduced failed on {} with a flat target, {} without",
                                   n, full_ok, red_ok, leg_ok, hold_ok, held_wells, dead_disagree, failed_grouped, failed_plain));
    for (const auto& [w, k] : dead_full) {
        BOOST_TEST_MESSAGE(fmt::format("  {} dead in full {} / reduced {} / legacy {} / reduced-hold {} of {}", w, k, dead_red[w], dead_leg[w], dead_hold[w], n));
    }
}

// The statement the whole design rests on, checked the other way round: an
// answer of the legacy rules, with the set it ended on, satisfies the full
// system's rows for that set. Node pressures, rates and bhps from the legacy
// answer go into the state; the controls it ended on go into the system;
// the residual must be at the legacy's own tolerance. The flattened group
// target is left out -- the legacy route has no group share.
BOOST_AUTO_TEST_CASE(legacy_answers_satisfy_the_equations)
{
    using Sys = DeckTrees::Sys;
    const char* dir = std::getenv("OPM_NETWORK_DUMP_PROD");
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    if (dir == nullptr || inc == nullptr || !std::filesystem::is_directory(dir)) {
        BOOST_TEST_MESSAGE("OPM_NETWORK_DUMP_PROD / OPM_VFP_INCLUDE not set, nothing to check");
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
    int n = 0, converged = 0, satisfied = 0, dead_ipr_rows = 0;
    double worst = 0.0;
    std::string worst_file;
    for (const auto& e : std::filesystem::directory_iterator(dir)) {
        if (e.path().extension() != ".txt") { continue; }
        std::ifstream in(e.path());
        std::string head; std::getline(in, head);
        if (head != "production") { continue; }
        auto [sys, guess] = NetworkSolve::readProduction<double>(in, props, units);
        sys.setGroupTarget(0.0);
        sys.setAnalyticJacobian(true);
        sys.setComplementarity(false);
        ++n;
        const auto rl = NetworkSolve::solveLegacy(sys, guess);
        if (!rl.converged) { continue; }
        ++converged;
        // The legacy answer as a state of the full system, its set as the controls.
        auto x = sys.start(rl.node_pressure);
        for (int w = 0; w < sys.numWells(); ++w) {
            const auto& well = sys.wells()[w];
            const double q = rl.well_rate[w];
            const double bhp = (q > 0.0 && well.ipr_b[1] < 0.0) ? (q - well.ipr_a[1]) / well.ipr_b[1] : well.bhp_limit;
            x[sys.wellBhpIndex(w)] = bhp;
            for (int ph = 0; ph < Sys::NP; ++ph) {
                x[sys.wellRateIndex(w, ph)] = (q > 0.0) ? std::max(sys.ipr(well, ph, bhp), 0.0) : 0.0;
            }
            const char c = rl.controls[w];
            sys.setControl(w, c == 'T' ? Sys::Control::Thp : c == 'B' ? Sys::Control::Bhp
                              : c == 'O' ? Sys::Control::OilRate : Sys::Control::Shut);
        }
        for (int nd = sys.numNodes(); nd >= 1; --nd) {
            for (int ph = 0; ph < Sys::NP; ++ph) {
                double sum = sys.nodeSource(nd)[ph];
                for (int w = 0; w < sys.numWells(); ++w) {
                    if (sys.wells()[w].node != nd) { continue; }
                    sum += sys.wells()[w].efficiency * (x[sys.wellRateIndex(w, ph)] + (ph == 2 ? sys.wells()[w].lift_gas : 0.0));
                }
                for (int c = 1; c <= sys.numNodes(); ++c) {
                    if (sys.nodes()[c].parent == nd) { sum += sys.nodes()[c].efficiency * x[sys.qIdx(c, ph)]; }
                }
                x[sys.qIdx(nd, ph)] = sum;
            }
        }
        auto r = sys.residual(x);
        // A dead well is q = 0 on every phase; the linear IPR rows cannot say
        // that (at zero oil the other lines are not at zero), so the full
        // system has no state for a dead well that satisfies them. Those
        // rows are the full system's defect, not the answer's; left out and
        // counted.
        for (int w = 0; w < sys.numWells(); ++w) {
            if (rl.controls[w] == 'D' || rl.controls[w] == 'S') {
                for (int ph = 0; ph < Sys::NP; ++ph) {
                    if (std::abs(r[sys.wellRateIndex(w, ph)]) > 0.02) { ++dead_ipr_rows; }
                    r[sys.wellRateIndex(w, ph)] = 0.0;
                }
            }
        }
        double m = 0.0;
        int at = -1;
        for (int i = 0; i < static_cast<int>(r.size()); ++i) {
            if (std::abs(r[i]) > m) { m = std::abs(r[i]); at = i; }
        }
        if (m < 0.02) { ++satisfied; }
        if (m > worst) {
            worst = m;
            std::string rows;
            std::vector<int> idx(r.size());
            std::iota(idx.begin(), idx.end(), 0);
            std::partial_sort(idx.begin(), idx.begin() + 4, idx.end(), [&](int a, int b) { return std::abs(r[a]) > std::abs(r[b]); });
            for (int k = 0; k < 4; ++k) {
                const int i = idx[k];
                std::string what = "?";
                if (i < sys.numNodes()) { what = "node " + sys.nodes()[i + 1].name; }
                else if (i < 4 * sys.numNodes()) { what = fmt::format("branch {} ph {}", sys.nodes()[(i - sys.numNodes()) / 3 + 1].name, (i - sys.numNodes()) % 3); }
                else {
                    for (int w = 0; w < sys.numWells(); ++w) {
                        for (int ph = 0; ph < 3; ++ph) { if (i == sys.wellRateIndex(w, ph)) { what = fmt::format("well {} ipr ph {} [{}]", sys.wells()[w].name, ph, rl.controls[w]); } }
                        if (i == sys.wellBhpIndex(w)) { what = fmt::format("well {} control [{}]", sys.wells()[w].name, rl.controls[w]); }
                    }
                }
                rows += fmt::format(" {}={:.3g} ({})", i, r[i], what);
            }
            worst_file = e.path().filename().string() + ":" + rows;
        }
    }
    BOOST_TEST_MESSAGE(fmt::format("{} dumps, legacy converged on {}: {} satisfy the full system's rows at its set"
                                   " (max scaled residual < 0.02, dead wells' IPR rows aside: {} of those violated);"
                                   " worst {:.3g} in {}",
                                   n, converged, satisfied, dead_ipr_rows, worst, worst_file));
    BOOST_CHECK_EQUAL(satisfied, converged);
}

// How far each flowing well is from its cliff: the highest node pressure
// at which its tubing still has a stable crossing, against the pressure it
// sits at. Diagnostic; run with OPM_NETWORK_DUMP_PROD and OPM_VFP_INCLUDE.
BOOST_AUTO_TEST_CASE(crossing_margin_on_the_dumps)
{
    const char* dir = std::getenv("OPM_NETWORK_DUMP_PROD");
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    if (dir == nullptr || inc == nullptr || !std::filesystem::is_directory(dir)) { return; }
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
    std::vector<std::filesystem::path> files;
    for (const auto& e : std::filesystem::directory_iterator(dir)) { if (e.path().extension() == ".txt") { files.push_back(e.path()); } }
    std::sort(files.begin(), files.end());
    for (const auto& f : files) {
        std::ifstream in(f);
        std::string head; std::getline(in, head);
        if (head != "production") { continue; }
        auto [sys, guess] = NetworkSolve::readProduction<double>(in, props, units);
        sys.setExactPotential(true);
        std::string line = f.filename().string() + ":";
        for (int w = 0; w < sys.numWells(); ++w) {
            const auto& well = sys.wells()[w];
            if (well.vfp_table <= 0) { continue; }
            const double p = guess[well.node];
            const double q = sys.thpPotential(well, p);
            // the cliff: bisect for the highest pressure with a crossing
            double lo = p, hi = p + 20.0 * unit::barsa;
            if (!(q > 0.0)) { hi = p; lo = std::max(p - 20.0 * unit::barsa, unit::barsa); }
            for (int it = 0; it < 30; ++it) {
                const double mid = 0.5 * (lo + hi);
                (sys.thpPotential(well, mid) > 0.0 ? lo : hi) = mid;
            }
            line += fmt::format("  {} p {:.2f} q {:.0f} start {:.0f} cliff {:.2f}{}", well.name, p / unit::barsa, q * 86400.0,
                                well.q_start * 86400.0, lo / unit::barsa, well.shut ? " (shut)" : "");
        }
        BOOST_TEST_MESSAGE(line);
    }
}

// The well solve as the well function. Each dumped well gets a curved inflow
// with the dump's linear IPR as its tangent at the operating point and a
// slope kappa*100 % different at shut-in; the reduced form re-linearises at
// every evaluation (solveWells). Against the fixed-IPR answer: how far the
// node pressures move, which wells change state, and what it costs.
BOOST_AUTO_TEST_CASE(the_well_solve_as_the_well_function)
{
    using Sys = DeckTrees::Sys;
    const char* dir = std::getenv("OPM_NETWORK_DUMP_PROD");
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    if (dir == nullptr || inc == nullptr || !std::filesystem::is_directory(dir)) { return; }
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
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    auto letters = [](const Sys& sys) { std::string s; for (int w = 0; w < sys.numWells(); ++w) { s += sys.controlLetter(w); } return s; };
    auto curved = [&](Sys& sys, const double kappa) {
        for (int w = 0; w < sys.numWells(); ++w) {
            const auto& well = sys.wells()[w];
            if (well.shut || well.vfp_table <= 0 || !(well.ipr_b[1] < 0.0)) { continue; }
            const double bhp0 = (well.q_start > 0.0) ? (well.q_start - well.ipr_a[1]) / well.ipr_b[1] : well.bhp_limit;
            const auto a = well.ipr_a, b = well.ipr_b;
            std::array<double, Sys::NP> c{};
            for (int ph = 0; ph < Sys::NP; ++ph) {
                const double shut_in = (b[ph] < 0.0) ? -a[ph] / b[ph] : bhp0 + 100.0 * unit::barsa;
                c[ph] = (shut_in > bhp0) ? kappa * b[ph] / (2.0 * (shut_in - bhp0)) : 0.0;
            }
            sys.setWellInflow(w, [a, b, c, bhp0](const double bhp) {
                std::array<double, Sys::NP> q{};
                for (int ph = 0; ph < Sys::NP; ++ph) { q[ph] = a[ph] + b[ph] * bhp + c[ph] * (bhp - bhp0) * (bhp - bhp0); }
                return q;
            });
        }
    };
    using clock = std::chrono::steady_clock;
    struct Tally { int n = 0, converged = 0, verified = 0, set_changed = 0, shut_changed = 0; double max_dp = 0.0, sum_dp = 0.0, ms = 0.0; long evals = 0, solves = 0; std::map<std::string, int> why; };
    std::map<std::string, Tally> tally;
    std::vector<std::filesystem::path> files;
    for (const auto& e : std::filesystem::directory_iterator(dir)) { if (e.path().extension() == ".txt") { files.push_back(e.path()); } }
    std::sort(files.begin(), files.end());
    int shown = 0;
    for (const auto& f : files) {
        std::ifstream in(f);
        std::string head; std::getline(in, head);
        if (head != "production") { continue; }
        auto [dumped, guess] = NetworkSolve::readProduction<double>(in, props, units);
        auto base = dumped; base.setAnalyticJacobian(true);
        const auto rb = NetworkSolve::solveReduced(base, guess, params, true);
        if (!rb.converged) { continue; }
        const std::string set_b = letters(base);
        for (const double kappa : {0.0, -0.3, 0.3}) {
            const std::string key = fmt::format("kappa {:+.1f}", kappa);
            auto& t = tally[key];
            auto sys = dumped; sys.setAnalyticJacobian(true);
            curved(sys, kappa);
            const auto t0 = clock::now();
            const auto r = NetworkSolve::solveReduced(sys, guess, params, true);
            t.ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
            ++t.n; t.evals += r.evaluations; t.solves += sys.wellSolves();
            if (!r.converged) { continue; }
            ++t.converged;
            const auto v = verifyAnswer(sys, r.node_pressure, r.well_rate, letters(sys));
            if (v.ok) { ++t.verified; } else { for (const auto& [k, cnt] : v.violations) { t.why[k] += cnt; } }
            double dp = 0.0;
            for (std::size_t i = 0; i < r.node_pressure.size(); ++i) { dp = std::max(dp, std::abs(r.node_pressure[i] - rb.node_pressure[i])); }
            t.max_dp = std::max(t.max_dp, dp); t.sum_dp += dp;
            const std::string set = letters(sys);
            if (set != set_b) { ++t.set_changed; }
            bool shut_diff = false;
            for (std::size_t i = 0; i < set.size(); ++i) { shut_diff |= (set[i] == 'S') != (set_b[i] == 'S'); }
            if (shut_diff) {
                ++t.shut_changed;
                if (shown++ < 6) {
                    BOOST_TEST_MESSAGE(fmt::format("  {} {}: set {} -> {}, node shift {:.3f} bar", f.filename().string(), key, set_b, set, dp / unit::barsa));
                }
            }
            if (kappa == 0.0) { BOOST_CHECK_SMALL(dp / unit::barsa, 1e-3); }
        }
    }
    for (const auto& [key, t] : tally) {
        BOOST_TEST_MESSAGE(fmt::format("  {}: of {}: converged {}, verified {}, set changed {}, shut set changed {}, node shift mean {:.3f} max {:.3f} bar, evals/solve {:.1f}, well solves/solve {:.0f}, ms/solve {:.2f}",
                                       key, t.n, t.converged, t.verified, t.set_changed, t.shut_changed,
                                       t.converged ? t.sum_dp / t.converged / unit::barsa : 0.0, t.max_dp / unit::barsa,
                                       t.n ? double(t.evals) / t.n : 0.0, t.n ? double(t.solves) / t.n : 0.0, t.n ? t.ms / t.n : 0.0));
        std::string why; for (const auto& [k, cnt] : t.why) { why += fmt::format(" {} x{}", k, cnt); }
        if (!why.empty()) { BOOST_TEST_MESSAGE("    not verified:" + why); }
    }
}

BOOST_AUTO_TEST_CASE(methods_on_the_dumps_scored)
{
    using Sys = DeckTrees::Sys;
    const char* dir = std::getenv("OPM_NETWORK_DUMP_PROD");
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    if (dir == nullptr || inc == nullptr || !std::filesystem::is_directory(dir)) {
        BOOST_TEST_MESSAGE("OPM_NETWORK_DUMP_PROD / OPM_VFP_INCLUDE not set, nothing to score");
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
    struct Score { int n = 0, converged = 0, verified = 0, wrong = 0, hysteresis = 0, shut_operable = 0, shut_checked = 0;
                   long lookups = 0; double ms = 0.0; std::map<std::string, int> why; };
    std::map<std::string, Score> scores;
    std::map<std::string, Score> scores_grouped;     // the dumps with a flat target
    int shown_wrong = 0;
    auto record = [&](const std::string& method, const bool grouped, const bool converged, const Verdict& v, const long lookups,
                      const double ms, const std::string& file = {}) {
        for (auto* sc : {&scores[method], grouped ? &scores_grouped[method] : nullptr}) {
            if (!sc) { continue; }
            ++sc->n; sc->lookups += lookups; sc->ms += ms;
            if (!converged) { continue; }
            ++sc->converged;
            if (v.hysteresis > 0) { ++sc->hysteresis; }
            if (v.ok) { ++sc->verified; } else { ++sc->wrong; for (const auto& [k, c] : v.violations) { sc->why[k] += c; } }
        }
        if (converged && !v.ok && method == "reduced" && shown_wrong++ < 6) {
            std::string why; for (const auto& [k, c] : v.violations) { why += fmt::format(" {} x{}", k, c); }
            BOOST_TEST_MESSAGE(fmt::format("  reduced wrong on {}:{}", file, why));
        }
    };
    using clock = std::chrono::steady_clock;
    auto ms_since = [](const clock::time_point t) { return std::chrono::duration<double, std::milli>(clock::now() - t).count(); };
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    auto letters = [](const Sys& sys) {
        std::string s;
        for (int w = 0; w < sys.numWells(); ++w) { s += sys.controlLetter(w); }
        return s;
    };
    for (const auto& e : std::filesystem::directory_iterator(dir)) {
        if (e.path().extension() != ".txt") { continue; }
        std::ifstream in(e.path());
        std::string head; std::getline(in, head);
        if (head != "production") { continue; }
        auto [dumped, guess] = NetworkSolve::readProduction<double>(in, props, units);
        const bool grouped = dumped.groupTarget() > 0.0;
        dumped.setAnalyticJacobian(true);
        const std::string fname = e.path().filename().string();
        {   // the legacy rules
            auto sys = dumped; sys.resetLookups();
            const auto t0 = clock::now();
            const auto r = NetworkSolve::solveLegacy(sys, guess);
            const double ms = ms_since(t0);
            const std::string c = r.controls;
            record("legacy rules", grouped, r.converged, r.converged ? verifyAnswer(sys, r.node_pressure, r.well_rate, c) : Verdict{}, sys.lookups(), ms);
        }
        {   // the legacy pressure update with the tree walk as its well model:
            // the simulator's fixed point on a system that has a group tree,
            // so the group rows can be judged against it too
            auto sys = dumped; sys.resetLookups();
            const auto t0 = clock::now();
            const auto r = NetworkSolve::solveLegacyOnWalk(sys, guess);
            const double ms = ms_since(t0);
            record("legacy rules + the walk", grouped, r.converged,
                   r.converged ? verifyAnswer(sys, r.node_pressure, r.well_rate, letters(sys)) : Verdict{},
                   sys.lookups(), ms);
        }
        {   // the reduced form, die rule -- and every shut well asked whether it could have flowed
            auto sys = dumped; sys.resetLookups();
            const auto t0 = clock::now();
            const auto r = NetworkSolve::solveReduced(sys, guess, params, true);
            const double ms = ms_since(t0);
            const auto v = r.converged ? verifyAnswer(sys, r.node_pressure, r.well_rate, letters(sys)) : Verdict{};
            record("reduced", grouped, r.converged, v, sys.lookups(), ms, fname);
            scores["reduced"].shut_checked += r.revived;   // reopened inside the solve, counted with the checks below
            if (r.converged && v.hysteresis > 0) {
                for (int w = 0; w < sys.numWells(); ++w) {
                    const auto& well = sys.wells()[w];
                    if (well.shut || sys.control(w) != Sys::Control::Shut || well.vfp_table <= 0) { continue; }
                    const double cross = sys.thpPotential(well, r.node_pressure[well.node]);
                    if (!(cross > 0.0)) { continue; }
                    // Give it its crossing as the flowing start and solve again
                    // from the settled pressures: is there an answer with it open?
                    auto again = dumped;
                    auto wells = again.wells();
                    Sys::Well opened = wells[w]; opened.q_start = cross;
                    // rebuild with the one well's start changed
                    Sys rebuilt(props, units);
                    for (int nd = 0; nd <= again.numNodes(); ++nd) { rebuilt.addNode(again.nodes()[nd], again.branchAlq(nd)); }
                    rebuilt.setTerminalPressure(again.terminalPressure());
                    for (int i = 0; i < again.numWells(); ++i) { rebuilt.addWell(i == w ? opened : wells[i]); }
                    rebuilt.setGroupTarget(again.groupTarget());
                    rebuilt.setAnalyticJacobian(true);
                    rebuilt.finish();
                    ++scores["reduced"].shut_checked;
                    const auto ra = NetworkSolve::solveReduced(rebuilt, r.node_pressure, params, true);
                    if (ra.converged && ra.well_rate[w] > 1e-9) {
                        const auto va = verifyAnswer(rebuilt, ra.node_pressure, ra.well_rate, letters(rebuilt));
                        if (va.ok && va.hysteresis == 0) {
                            ++scores["reduced"].shut_operable;
                            if (scores["reduced"].shut_operable <= 4) {
                                BOOST_TEST_MESSAGE(fmt::format("  {} shut but operable: {} flows at {:.0f} sm3/d, node {:.2f} -> {:.2f} bar",
                                                               fname, well.name, ra.well_rate[w] * 86400.0,
                                                               convert::to(r.node_pressure[well.node], bars),
                                                               convert::to(ra.node_pressure[well.node], bars)));
                            }
                        }
                    }
                }
            }
        }
        for (const bool cmpl : {false, true}) {
            for (const bool dead_rule : {false, true}) {
                auto sys = dumped; sys.resetLookups(); sys.setComplementarity(cmpl); sys.setDeadWhenCannotLift(dead_rule);
                if (dead_rule) { sys.resetDead(); }
                const auto t0 = clock::now();
                const auto r = NetworkSolve::solve(sys, guess, params, NetworkSolve::FullStep{});
                const double ms = ms_since(t0);
                const std::string name = std::string(cmpl ? "full complementarity" : "full active set") + (dead_rule ? " + dead rule" : "");
                record(name, grouped, r.converged, r.converged ? verifyAnswer(sys, r.node_pressure, r.well_rate, letters(sys)) : Verdict{}, sys.lookups(), ms);
            }
        }
    }
    for (const auto* tab : {&scores, &scores_grouped}) {
        const std::string heading = (tab == &scores) ? "all dumps:" : "dumps with a flat group target:";
        BOOST_TEST_MESSAGE(heading);
        for (const auto& [m, sc] : *tab) {
            std::string why;
            std::vector<std::pair<int, std::string>> top;
            for (const auto& [k, c] : sc.why) { top.emplace_back(c, k); }
            std::sort(top.rbegin(), top.rend());
            for (std::size_t i = 0; i < top.size() && i < 3; ++i) { why += fmt::format(" {} x{}", top[i].second, top[i].first); }
            BOOST_TEST_MESSAGE(fmt::format("  {:34} of {}: converged {:4}, verified {:4} ({:3} shut-stays-shut{}), wrong {:4}, lookups/solve {:6.0f}, ms/solve {:5.2f};{}",
                                           m, sc.n, sc.converged, sc.verified, sc.hysteresis,
                                           sc.shut_checked ? fmt::format(", {} of {} shut wells operable", sc.shut_operable, sc.shut_checked) : std::string{},
                                           sc.wrong, double(sc.lookups) / std::max(sc.n, 1), sc.ms / std::max(sc.n, 1), why));
        }
    }
}

// The routes on layered instances: deeper networks and trees, targets on
// several levels at once, and a quarter of the wells weak enough that the
// tubing caps them below their share so the tree hands the difference to
// the others. Every answer through the same judge.
BOOST_AUTO_TEST_CASE(generated_layered_cases_scored)
{
    using Sys = DeckTrees::Sys;
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    struct Score { int n = 0, converged = 0, verified = 0, wrong = 0, oracle_agrees = 0, oracle_off = 0, oracle_violates = 0; long it = 0, lookups = 0; double ms = 0.0; std::map<std::string, int> why; };
    std::map<std::string, Score> scores;
    using clock = std::chrono::steady_clock;
    auto record = [&](const std::string& m, const bool converged, const Verdict& v, const int it, const long lookups, const double ms) {
        auto& sc = scores[m];
        ++sc.n; sc.it += it; sc.lookups += lookups; sc.ms += ms;
        if (!converged) { return; }
        ++sc.converged;
        if (v.ok) { ++sc.verified; } else { ++sc.wrong; for (const auto& [k, c] : v.violations) { sc.why[k] += c; } }
    };
    // Stein's balancer, fed the capacities at an answer's pressures, against
    // that answer's rates.
    std::string instance;
    auto oracle = [&](const std::string& m, Sys& sys, const std::vector<double>& state, const std::vector<double>& rates, const Schedule& sched) {
        const auto stein = steinAllocation(sys, state, sched, 0);
        if (stein.empty()) { return; }
        double worst = 0.0;
        for (int w = 0; w < sys.numWells(); ++w) {
            worst = std::max(worst, std::abs(rates[w] * 86400.0 - stein[w]) / std::max(std::abs(stein[w]), 1.0));
        }
        if (worst < 0.005) { ++scores[m].oracle_agrees; return; }
        // His allocation above some group's own limit is the transparent-group
        // walk-through reported to him, not a disagreement about the answer.
        const auto& groups = sys.groups();
        bool violates = false;
        for (int g = 0; g < sys.numGroups() && !violates; ++g) {
            if (!(groups[g].target > 0.0)) { continue; }
            const auto c = Sys::modeWeights(groups[g].mode, {});
            double on = 0.0;
            for (int w = 0; w < sys.numWells(); ++w) {
                bool under = false;
                for (int a = sys.wells()[w].group; a >= 0; a = groups[a].parent) { if (a == g) { under = true; break; } }
                if (!under) { continue; }
                const auto& well = sys.wells()[w];
                const double q_oil = stein[w] / 86400.0;
                const double bhp = (q_oil - well.ipr_a[1]) / well.ipr_b[1];
                for (int ph = 0; ph < Sys::NP; ++ph) { on += c[ph] * std::max(well.ipr_a[ph] + well.ipr_b[ph] * bhp, 0.0); }
            }
            if (on > groups[g].target * 1.005) {
                violates = true;
                if (m == "C reduced") {
                    BOOST_TEST_MESSAGE(fmt::format("    {}: Stein puts {:.0f} on {} against its {:.0f} (has guide rate: {})",
                                                   instance, on * 86400.0, groups[g].name, groups[g].target * 86400.0,
                                                   sched.getGroup(groups[g].name, 0).productionControls({}).guide_rate > 0.0 ? "yes" : "no"));
                }
            }
        }
        if (!violates && m == "C reduced") {
            int wm = 0; double dm = 0.0;
            for (int w = 0; w < sys.numWells(); ++w) {
                const double d = std::abs(rates[w] * 86400.0 - stein[w]) / std::max(std::abs(stein[w]), 1.0);
                if (d > dm) { dm = d; wm = w; }
            }
            BOOST_TEST_MESSAGE(fmt::format("    {}: Stein differs, worst {}: ours {:.1f}, his {:.1f} sm3/d (control {})",
                                           instance, sys.wells()[wm].name, rates[wm] * 86400.0, stein[wm], sys.controlLetter(wm)));
        }
        (violates ? scores[m].oracle_violates : scores[m].oracle_off) += 1;
    };
    auto letters = [](const Sys& sys) { std::string s; for (int w = 0; w < sys.numWells(); ++w) { s += sys.controlLetter(w); } return s; };
    struct Shape { int wells, nodes, groups, gdepth, ndepth, seeds; };
    const std::vector<Shape> shapes{{30, 12, 10, 5, 5, 4}, {80, 30, 20, 6, 5, 3}, {160, 50, 30, 6, 6, 2},
                                    {120, 60, 60, 10, 8, 2}, {300, 80, 80, 12, 10, 1}};
    int shown = 0;
    const char* only = std::getenv("OPM_LAYERED_ONLY");     // "shape,seed" to run one instance
    for (std::size_t si = 0; si < shapes.size(); ++si) {
        const auto& sh = shapes[si];
        for (int seed = 1; seed <= sh.seeds; ++seed) {
            if (only && fmt::format("{},{}", si, seed) != only) { continue; }
            GenSpec spec{sh.wells, sh.nodes, sh.groups, static_cast<unsigned>(seed)};
            spec.group_depth = sh.gdepth; spec.net_depth = sh.ndepth; spec.stiff = 0.25; spec.target_fraction = 0.5;
            std::string what;
            const auto text = generateDeck(spec, &what);
            instance = fmt::format("shape {} seed {}", si, seed);
            if (const char* dump_dir = std::getenv("OPM_LAYERED_DUMP")) {   // the deck text, to reproduce elsewhere
                std::ofstream(std::filesystem::path(dump_dir) / fmt::format("layered_shape{}_seed{}.DATA", si, seed)) << text;
            }
            DeckTrees dt(text, DeckTrees::FromText{});
            DeckTrees::Ipr ipr; ipr.j_scale = 2.5; ipr.j_of = weakWells(spec);
            // node_order, and so the guess, exist only after a build.
            const auto guess = [&] { (void)dt.build(0, ipr); return dt.guess(20.0); }();
            {   // B2, the set frozen between converged solves
                auto sys = dt.build(0, ipr); sys.resetLookups();
                const auto t0 = clock::now();
                const auto r = NetworkSolve::solveWithTree(sys, guess, params, NetworkSolve::FullStep{});
                const double ms = std::chrono::duration<double, std::milli>(clock::now() - t0).count();
                const auto v = r.result.converged ? verifyAnswer(sys, r.result.node_pressure, r.result.well_rate, letters(sys)) : Verdict{};
                record("B2 outer loop", r.result.converged && r.consistent, v, r.inner_iterations, sys.lookups(), ms);
                if (r.result.converged) { oracle("B2 outer loop", sys, r.result.state, r.result.well_rate, *dt.schedule); }
                if (r.result.converged && !v.ok && shown++ < 4) {
                    std::string why; for (const auto& [k, c] : v.violations) { why += fmt::format(" {} x{}", k, c); }
                    BOOST_TEST_MESSAGE(fmt::format("  B2 on {} seed {}: {}", what, seed, why));
                }
            }
            {   // C, the reduced form
                auto sys = dt.build(0, ipr); sys.resetLookups();
                const auto t0 = clock::now();
                const auto r = NetworkSolve::solveReduced(sys, guess, params, true);
                const double ms = std::chrono::duration<double, std::milli>(clock::now() - t0).count();
                record("C reduced", r.converged, r.converged ? verifyAnswer(sys, r.node_pressure, r.well_rate, letters(sys)) : Verdict{}, r.iterations, sys.lookups(), ms);
                if (r.converged) { oracle("C reduced", sys, sys.reducedState(), r.well_rate, *dt.schedule); }
            }
            {   // B1, the set every iteration
                auto sys = dt.build(0, ipr); sys.resetLookups();
                sys.setGroupActiveSet(true);
                const auto t0 = clock::now();
                const auto r = NetworkSolve::solve(sys, guess, params, NetworkSolve::FullStep{});
                const double ms = std::chrono::duration<double, std::milli>(clock::now() - t0).count();
                record("B1 every iteration", r.converged, r.converged ? verifyAnswer(sys, r.node_pressure, r.well_rate, letters(sys)) : Verdict{}, r.iterations, sys.lookups(), ms);
            }
            {   // A, Fischer-Burmeister
                auto sys = dt.build(0, ipr); sys.resetLookups();
                const auto t0 = clock::now();
                const auto r = NetworkSolve::solve(sys, guess, params, NetworkSolve::FullStep{});
                const double ms = std::chrono::duration<double, std::milli>(clock::now() - t0).count();
                record("A complementarity", r.converged, r.converged ? verifyAnswer(sys, r.node_pressure, r.well_rate, letters(sys)) : Verdict{}, r.iterations, sys.lookups(), ms);
            }
            {   // the legacy fixed point, with the walk as its well model: the
                // simulator's damped, capped pressure update on the reduced residual
                auto sys = dt.build(0, ipr); sys.resetLookups();
                const auto t0 = clock::now();
                NetworkSolve::LegacyParameters<double> lp; lp.max_iterations = 3000;
                const auto r = NetworkSolve::solveLegacyOnWalk(sys, guess, lp);
                const double ms = std::chrono::duration<double, std::milli>(clock::now() - t0).count();
                record("legacy fixed point on the walk", r.converged, r.converged ? verifyAnswer(sys, r.node_pressure, r.well_rate, letters(sys)) : Verdict{}, r.iterations, sys.lookups(), ms);
            }
            BOOST_TEST_MESSAGE(fmt::format("seed {}: {}", seed, what));
        }
    }
    for (const auto& [m, sc] : scores) {
        std::string why;
        std::vector<std::pair<int, std::string>> top;
        for (const auto& [k, c] : sc.why) { top.emplace_back(c, k); }
        std::sort(top.rbegin(), top.rend());
        for (std::size_t i = 0; i < top.size() && i < 3; ++i) { why += fmt::format(" {} x{}", top[i].second, top[i].first); }
        BOOST_TEST_MESSAGE(fmt::format("  {:30} of {}: converged {:2}, verified {:2}, wrong {:2}, Stein agrees/off/above-a-limit {}/{}/{}, it/instance {:5.1f}, lookups {:7.0f}, ms {:7.1f};{}",
                                       m, sc.n, sc.converged, sc.verified, sc.wrong, sc.oracle_agrees, sc.oracle_off, sc.oracle_violates,
                                       double(sc.it) / std::max(sc.n, 1), double(sc.lookups) / std::max(sc.n, 1), sc.ms / std::max(sc.n, 1), why));
    }
}

BOOST_AUTO_TEST_CASE(generated_case_anatomy)
{
    using Sys = DeckTrees::Sys;
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    int shown = 0;
    for (int seed = 1; seed <= 60 && shown < 3; ++seed) {
        std::string what;
        const auto text = generateDeck({20, 6, 6, static_cast<unsigned>(seed)}, &what);
        DeckTrees dt(text, DeckTrees::FromText{});
        DeckTrees::Ipr ipr; ipr.j_scale = 2.5;
        auto reduced = dt.build(0, ipr);
        auto outer = dt.build(0, ipr);
        const auto rr = NetworkSolve::solveReduced(reduced, dt.guess(20.0), params, true);
        const auto ro = NetworkSolve::solveWithTree(outer, dt.guess(20.0), params, NetworkSolve::FullStep{});
        std::vector<double> stein;
        double gap = 0.0;
        if (rr.converged) {
            stein = steinAllocation(reduced, reduced.reducedState(), *dt.schedule, 0);
            for (int w = 0; w < reduced.numWells() && !stein.empty(); ++w) {
                gap = std::max(gap, std::abs(rr.well_rate[w] * 86400.0 - stein[w]) / std::max(std::abs(stein[w]), 1.0));
            }
        }
        bool off = false;
        if (rr.converged && ro.result.converged) {
            for (int w = 0; w < reduced.numWells(); ++w) {
                if (std::abs(ro.result.well_rate[w] - rr.well_rate[w]) > 0.005 * std::max(rr.well_rate[w], 1.0 / 86400.0)) { off = true; }
            }
        }
        if (rr.converged && gap <= 0.005 && !off && ro.result.converged) { continue; }
        ++shown;
        BOOST_TEST_MESSAGE(fmt::format("seed {}: {}; reduced {}{} in {} it, oracle gap {:.3g} %; outer {} {} passes{}{}",
                                       seed, what, rr.converged ? "converged" : "FAILED",
                                       stein.empty() ? " (oracle rejected)" : "", rr.iterations, 100 * gap,
                                       ro.result.converged ? "converged" : "FAILED", ro.passes,
                                       ro.consistent ? "" : " inconsistent", off ? "; DISAGREE" : ""));
        // The state the reduced form stopped at, with its set.
        reduced.updateControls(reduced.reducedState());
        const auto& groups = reduced.groups();
        std::string gl;
        for (int g = 0; g < reduced.numGroups(); ++g) {
            int held = 0;
            for (int c = 0; c < reduced.numGroups(); ++c) {
                if (groups[c].parent == g && reduced.groupBind(c) == Sys::GroupBind::Share) { ++held; }
            }
            for (int w = 0; w < reduced.numWells(); ++w) {
                if (reduced.wells()[w].group == g && reduced.control(w) == Sys::Control::Tree) { ++held; }
            }
            const char bind = reduced.groupBind(g) == Sys::GroupBind::Free ? 'F'
                            : reduced.groupBind(g) == Sys::GroupBind::Own ? 'O' : 'S';
            gl += fmt::format("\n    {:5} <- {:5} {} target {:.0f} guide {:.0f} bind {} held-children {}{}",
                              groups[g].name, groups[g].parent >= 0 ? groups[groups[g].parent].name : "-",
                              groups[g].mode == Sys::Mode::Liquid ? "LRAT" : "ORAT",
                              groups[g].target * 86400.0, groups[g].guide * 86400.0, bind, held,
                              (bind != 'F' && held == 0) ? "  <-- lambda without a column" : "");
        }
        BOOST_TEST_MESSAGE("  groups:" << gl);
        std::string nl;
        for (std::size_t n = 0; n < dt.node_order.size(); ++n) {
            nl += fmt::format(" {}({:.1f})", dt.node_order[n], convert::to(rr.node_pressure[n], bars));
        }
        BOOST_TEST_MESSAGE("  nodes:" << nl);
        std::string wl;
        for (int w = 0; w < reduced.numWells(); ++w) {
            const auto& well = reduced.wells()[w];
            wl += fmt::format("\n    {:4} in {:5} at {:5} cap {:7.1f} ctl {} reduced {:7.1f} outer {:7.1f} stein {:7.1f}",
                              well.name, groups[well.group].name, dt.node_order[well.node],
                              reduced.ownAllowance(w) * 86400.0, reduced.controlLetter(w),
                              rr.well_rate[w] * 86400.0,
                              ro.result.converged ? ro.result.well_rate[w] * 86400.0 : 0.0,
                              stein.empty() ? 0.0 : stein[w]);
        }
        BOOST_TEST_MESSAGE("  wells:" << wl);
    }
    BOOST_TEST_MESSAGE("shown " << shown << " of the first 60 seeds");
}

BOOST_AUTO_TEST_CASE(production_step_bounds)
{
    ProductionCase c;
    auto system = c.system();
    auto x = system.start(ProductionCase::guess());
    std::vector<double> dx(system.size(), 0.0);
    dx[system.pIdx(1)] = convert::from(1000.0, bars);
    auto out = system.limitStep(x, dx);
    BOOST_TEST_MESSAGE("1000 bar step -> " << convert::to(out[system.pIdx(1)], bars) << " bar");
    BOOST_CHECK_LE(out[system.pIdx(1)], convert::from(50.0, bars) * 1.0001);
    dx[system.pIdx(1)] = -x[system.pIdx(1)] - convert::from(5.0, bars);
    out = system.limitStep(x, dx);
    BOOST_CHECK_GE(x[system.pIdx(1)] + out[system.pIdx(1)], convert::from(1.0, bars) * 0.99);
}

// The tubing curve continued below its loading hump so that every inflow
// line crosses it once (Stein's suggestion): solve with no shut logic in the
// loop, then shut every well whose converged point is on the continuation,
// all at once, and solve again until none is. Judged on the real tables.
// Run with OPM_NETWORK_DUMP_PROD and OPM_VFP_INCLUDE; OPM_EXTENSION_TRACE
// prints the wells of the first systems that differ from the baseline.
BOOST_AUTO_TEST_CASE(the_tubing_extension_on_the_dumps)
{
    using Sys = DeckTrees::Sys;
    const char* dir = std::getenv("OPM_NETWORK_DUMP_PROD");
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    if (dir == nullptr || inc == nullptr || !std::filesystem::is_directory(dir)) { return; }
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
    static const bool trace = std::getenv("OPM_EXTENSION_TRACE") != nullptr;
    struct Score { int n = 0, converged = 0, verified = 0, wrong = 0, hysteresis = 0, closed = 0;
                   int passes1 = 0, passes2 = 0, passes3 = 0, same_shut = 0, more_shut = 0, fewer_shut = 0, other_shut = 0;
                   int rescued = 0, cliff_resolved = 0, compared = 0, on_cliff = 0, cliffs = 0, stalls = 0, revived = 0, set_changes = 0;
                   long iterations = 0, evaluations = 0, lookups = 0; double ms = 0.0, max_dp = 0.0, sum_dp = 0.0;
                   std::map<std::string, int> why; };
    std::map<std::string, Score> scores;
    std::array<int, 3> gap_bins{};   // closed with a gap < 0.1, < 1, >= 1 bar
    // The closing decision: three strategies against closing all at once,
    // and a re-open band measured on the wells the rule closed.
    struct ClosingScore { int n = 0, converged = 0, passes1 = 0, passes2 = 0, passes3 = 0, differs = 0, more = 0, fewer = 0; };
    std::map<std::string, ClosingScore> closing_scores;
    std::array<int, 4> dip_bins{};              // shut wells: line dips below the curve by < 0.1 / < 0.5 / >= 0.5 bar / not at all
    const std::array<double, 3> bands{0.0, 0.1, 0.5};
    std::array<int, 3> reopen_candidates{}, reopen_stable{};
    using clock = std::chrono::steady_clock;
    auto ms_since = [](const clock::time_point t) { return std::chrono::duration<double, std::milli>(clock::now() - t).count(); };
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    auto letters = [](const Sys& sys) { std::string s; for (int w = 0; w < sys.numWells(); ++w) { s += sys.controlLetter(w); } return s; };
    std::vector<std::filesystem::path> files;
    for (const auto& e : std::filesystem::directory_iterator(dir)) { if (e.path().extension() == ".txt") { files.push_back(e.path()); } }
    std::sort(files.begin(), files.end());
    int shown = 0, traced = 0;
    for (const auto& f : files) {
        std::ifstream in(f);
        std::string head; std::getline(in, head);
        if (head != "production") { continue; }
        auto [dumped, guess] = NetworkSolve::readProduction<double>(in, props, units);
        dumped.setAnalyticJacobian(true);
        const std::string fname = f.filename().string();
        // Baseline: the reduced form with the well model's dead rule.
        std::string set0; bool conv0 = false, ok0 = false, cliff0 = false; std::vector<double> p0;
        {
            auto sys = dumped; sys.resetLookups();
            const auto t0 = clock::now();
            const auto r = NetworkSolve::solveReduced(sys, guess, params, true);
            auto& sc = scores["reduced (baseline)"];
            ++sc.n; sc.ms += ms_since(t0); sc.lookups += sys.lookups(); sc.iterations += r.iterations; sc.evaluations += r.evaluations; ++sc.passes1;
            sc.cliffs += r.cliffs; sc.stalls += r.stalls; sc.revived += r.revived; sc.set_changes += r.set_changes;
            conv0 = r.converged; cliff0 = r.on_cliff; set0 = letters(sys); p0 = r.node_pressure;
            if (!conv0 && trace) { BOOST_TEST_MESSAGE(fmt::format("  baseline did not converge on {}: residual {:.3g} after {} it, set {}", fname, r.residual, r.iterations, set0)); }
            if (r.on_cliff) { ++sc.on_cliff; }
            if (conv0) {
                ++sc.converged;
                const auto v = verifyAnswer(sys, r.node_pressure, r.well_rate, set0);
                ok0 = v.ok;
                if (v.hysteresis > 0) { ++sc.hysteresis; }
                if (v.ok) { ++sc.verified; } else { ++sc.wrong; for (const auto& [k, c] : v.violations) { sc.why[k] += c; } }
            }
        }
        {
            auto sys = dumped; sys.resetLookups();
            const auto t0 = clock::now();
            const auto rx = NetworkSolve::solveReducedOnExtension(sys, guess, params);
            auto& sc = scores["extension"];
            ++sc.n; sc.ms += ms_since(t0); sc.lookups += sys.lookups(); sc.iterations += rx.iterations; sc.evaluations += rx.evaluations; sc.closed += rx.closed;
            (rx.passes == 1 ? sc.passes1 : rx.passes == 2 ? sc.passes2 : sc.passes3) += 1;
            sc.cliffs += rx.last.cliffs; sc.stalls += rx.last.stalls; sc.revived += rx.last.revived; sc.set_changes += rx.last.set_changes;
            if (!rx.converged && trace) { BOOST_TEST_MESSAGE(fmt::format("  extension did not converge on {}: residual {:.3g} after {} it, pass {}, set {}, baseline {}{}", fname, rx.last.residual, rx.last.iterations, rx.passes, letters(sys), set0, conv0 ? "" : " (baseline neither)")); }
            for (std::size_t w = 0; w < rx.closed_wells.size(); ++w) {
                if (!rx.closed_wells[w]) { continue; }
                const double g = convert::to(rx.closed_gap[w], bars);
                ++gap_bins[g < 0.1 ? 0 : g < 1.0 ? 1 : 2];
            }
            if (!rx.converged) { continue; }
            ++sc.converged;
            if (rx.last.on_cliff) { ++sc.on_cliff; }
            const std::string set = letters(sys);
            const auto v = verifyAnswer(sys, rx.last.node_pressure, rx.last.well_rate, set);
            if (v.hysteresis > 0) { ++sc.hysteresis; }
            if (v.ok) { ++sc.verified; } else { ++sc.wrong; for (const auto& [k, c] : v.violations) { sc.why[k] += c; } }
            if (trace && (!v.ok || set != set0) && traced++ < 5) {
                std::string why; for (const auto& [k, c] : v.violations) { why += fmt::format(" {} x{}", k, c); }
                BOOST_TEST_MESSAGE(fmt::format("  trace {}: baseline {} extension {} passes {} closed {}{}", fname, set0, set, rx.passes, rx.closed, why));
                for (int w = 0; w < sys.numWells(); ++w) {
                    const auto& well = sys.wells()[w];
                    if (well.vfp_table <= 0) { continue; }
                    const double pn = well.node == 0 ? sys.terminalPressure() : rx.last.node_pressure[well.node];
                    const double pg = well.node == 0 ? sys.terminalPressure() : guess[well.node];
                    const double shut_in = well.ipr_b[1] < 0.0 ? -well.ipr_a[1] / well.ipr_b[1] : 0.0;
                    const double po = sys.thpPotential(well, pn), pog = sys.thpPotential(well, pg);
                    sys.setTubingExtension(true);
                    const double pe = sys.thpPotential(well, pn), peg = sys.thpPotential(well, pg);
                    sys.setTubingExtension(false);
                    auto show = [](const double q) { return q > 1e29 ? std::string("nobind") : fmt::format("{:.1f}", q * 86400.0); };
                    BOOST_TEST_MESSAGE(fmt::format("    {} {}->{}: q {:.1f} sm3/d at node {:.3f} bar (guess {:.3f}); crossing table {} continued {}; at guess {} / {}; q_start {:.1f} shut-in {:.1f} limit {:.1f} bar{}",
                                                   well.name, set0[w], set[w], rx.last.well_rate[w] * 86400.0, convert::to(pn, bars), convert::to(pg, bars),
                                                   show(po), show(pe), show(pog), show(peg), well.q_start * 86400.0, convert::to(shut_in, bars), convert::to(well.bhp_limit, bars),
                                                   rx.closed_wells[w] ? fmt::format(", closed with gap {:.3f} bar", convert::to(rx.closed_gap[w], bars)) : std::string{}));
                }
            }
            if (!conv0) { ++sc.rescued; }
            if (cliff0 && !rx.last.on_cliff) { ++sc.cliff_resolved; }
            if (conv0) {
                ++sc.compared;
                int more = 0, fewer = 0;
                for (std::size_t w = 0; w < set.size(); ++w) {
                    if (set[w] == 'S' && set0[w] != 'S') { ++more; }
                    if (set[w] != 'S' && set0[w] == 'S') { ++fewer; }
                }
                if (more == 0 && fewer == 0) { ++sc.same_shut; } else if (fewer == 0) { ++sc.more_shut; } else if (more == 0) { ++sc.fewer_shut; } else { ++sc.other_shut; }
                if (ok0 && v.ok) {
                    double dp = 0.0;
                    for (std::size_t i = 1; i < p0.size() && i < rx.last.node_pressure.size(); ++i) { dp = std::max(dp, std::abs(p0[i] - rx.last.node_pressure[i])); }
                    sc.max_dp = std::max(sc.max_dp, dp); sc.sum_dp += dp;
                }
                if ((more || fewer) && shown++ < 8) {
                    BOOST_TEST_MESSAGE(fmt::format("  {}: baseline {} -> {} ({} passes, {} closed by the rule{})",
                                                   fname, set0, set, rx.passes, rx.closed, v.ok ? "" : ", judged wrong"));
                }
            }
        }
        {
            using NetworkSolve::Closing;
            std::string set_all;
            for (const auto& [closing, name] : {std::pair{Closing::All, "all at once"}, std::pair{Closing::Sequential, "one at a time"}, std::pair{Closing::Tiered, "tiered"}}) {
                auto sys = dumped;
                const auto rx = NetworkSolve::solveReducedOnExtension(sys, guess, params, closing);
                auto& cs = closing_scores[name];
                ++cs.n;
                (rx.passes == 1 ? cs.passes1 : rx.passes == 2 ? cs.passes2 : cs.passes3) += 1;
                if (!rx.converged) { continue; }
                ++cs.converged;
                const std::string set = letters(sys);
                if (closing == Closing::All) {
                    set_all = set;
                    // The re-open band: how far below the curve each closed
                    // well's line dips at the settled pressure, and whether
                    // re-opening it gives an answer with it open.
                    for (int w = 0; w < sys.numWells(); ++w) {
                        if (!rx.closed_wells[w]) { continue; }
                        const auto& well = sys.wells()[w];
                        std::array<double, Sys::NP> q{};
                        for (int ph = 0; ph < Sys::NP; ++ph) { q[ph] = std::max(sys.ipr(well, ph, rx.closed_bhp[w]), 0.0); }
                        const double pn = well.node == 0 ? sys.terminalPressure() : rx.last.node_pressure[well.node];
                        const auto tp = sys.touchingPoint(well, pn, q);
                        const double dip = tp.valid ? convert::to(-tp.gap, bars) : -1.0;
                        ++dip_bins[dip <= 0.0 ? 3 : dip < 0.1 ? 0 : dip < 0.5 ? 1 : 2];
                        const double cross = sys.thpPotential(well, pn);
                        for (std::size_t b = 0; b < bands.size(); ++b) {
                            if (!(dip > bands[b]) || !(cross > 0.0) || cross > 1e29) { continue; }
                            ++reopen_candidates[b];
                            auto again = sys;
                            again.reviveWell(w, cross);
                            const auto r2 = NetworkSolve::solveReducedOnExtension(again, rx.last.node_pressure, params, Closing::All, 20, /*keep_dead=*/true);
                            if (r2.converged && letters(again)[w] != 'S') { ++reopen_stable[b]; }
                        }
                    }
                } else if (!set_all.empty()) {
                    int more = 0, fewer = 0;
                    for (std::size_t w = 0; w < set.size(); ++w) {
                        if (set[w] == 'S' && set_all[w] != 'S') { ++more; }
                        if (set[w] != 'S' && set_all[w] == 'S') { ++fewer; }
                    }
                    if (more || fewer) { ++cs.differs; cs.more += more; cs.fewer += fewer; }
                }
            }
        }
    }
    for (const auto& [m, cs] : closing_scores) {
        BOOST_TEST_MESSAGE(fmt::format("  closing {:14} of {}: converged {:4}, passes 1/2/3+: {}/{}/{}, shut set differs from all-at-once in {} systems (wells: {} more shut, {} fewer)",
                                       m, cs.n, cs.converged, cs.passes1, cs.passes2, cs.passes3, cs.differs, cs.more, cs.fewer));
    }
    BOOST_TEST_MESSAGE(fmt::format("  re-open band: of the wells the rule closed, at the settled pressure the line dips below the curve by < 0.1 / < 0.5 / >= 0.5 bar in {} / {} / {}, not at all in {}; candidates with a band of 0 / 0.1 / 0.5 bar: {} / {} / {}, staying open once re-opened: {} / {} / {}",
                                   dip_bins[0], dip_bins[1], dip_bins[2], dip_bins[3],
                                   reopen_candidates[0], reopen_candidates[1], reopen_candidates[2], reopen_stable[0], reopen_stable[1], reopen_stable[2]));
    for (const auto& [m, sc] : scores) {
        std::string why;
        std::vector<std::pair<int, std::string>> top;
        for (const auto& [k, c] : sc.why) { top.emplace_back(c, k); }
        std::sort(top.rbegin(), top.rend());
        for (std::size_t i = 0; i < top.size() && i < 3; ++i) { why += fmt::format(" {} x{}", top[i].second, top[i].first); }
        BOOST_TEST_MESSAGE(fmt::format("  {:20} of {}: converged {:4}, verified {:4} ({:3} shut-stays-shut), wrong {:4}, on a cliff {:3}, cliff cutbacks {:3}, stalls {:3}, revived {:3}, set changes {:4}, it/solve {:4.1f}, evals/solve {:5.1f}, lookups/solve {:6.0f}, ms/solve {:5.2f};{}",
                                       m, sc.n, sc.converged, sc.verified, sc.hysteresis, sc.wrong, sc.on_cliff, sc.cliffs, sc.stalls, sc.revived, sc.set_changes,
                                       double(sc.iterations) / std::max(sc.n, 1), double(sc.evaluations) / std::max(sc.n, 1), double(sc.lookups) / std::max(sc.n, 1), sc.ms / std::max(sc.n, 1), why));
        if (m == "reduced (baseline)") { continue; }
        BOOST_TEST_MESSAGE(fmt::format("  {:20} closed wells, line missed the curve by < 0.1 / < 1 / >= 1 bar: {} / {} / {}", "", gap_bins[0], gap_bins[1], gap_bins[2]));
        BOOST_TEST_MESSAGE(fmt::format("  {:20} passes 1/2/3+: {}/{}/{}, wells closed by the rule {}; vs baseline ({} compared): same shut set {}, more shut {}, fewer shut {}, mixed {}; rescued {}, cliffs resolved {}; node pressure max {:.3f} mean {:.3f} bar",
                                       "", sc.passes1, sc.passes2, sc.passes3, sc.closed, sc.compared, sc.same_shut, sc.more_shut, sc.fewer_shut, sc.other_shut,
                                       sc.rescued, sc.cliff_resolved, convert::to(sc.max_dp, bars), convert::to(sc.sum_dp / std::max(sc.compared, 1), bars)));
    }
}

// The IPR of a well that cannot lift: where its tangent is taken decides
// whether the linear rows say it lifts. Each dumped well gets a curved
// inflow (the dump's line as its tangent at the operating point, kappa*100 %
// off in slope at shut-in), and the shut decision of the solve is checked
// against the exact one for the curve at the settled pressure: the curve
// misses the tubing everywhere, or not. Tangent rule off (where the well
// last flowed, or its bhp limit) against on (the touching point, or a rate
// limit below the crossing), with and without the tubing continuation.
// Run with OPM_NETWORK_DUMP_PROD and OPM_VFP_INCLUDE.
BOOST_AUTO_TEST_CASE(the_shut_well_ipr_on_the_dumps)
{
    using Sys = DeckTrees::Sys;
    const char* dir = std::getenv("OPM_NETWORK_DUMP_PROD");
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    if (dir == nullptr || inc == nullptr || !std::filesystem::is_directory(dir)) { return; }
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
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    static const bool trace = std::getenv("OPM_EXTENSION_TRACE") != nullptr;
    auto letters = [](const Sys& sys) { std::string s; for (int w = 0; w < sys.numWells(); ++w) { s += sys.controlLetter(w); } return s; };
    using Inflow = std::function<std::array<double, Sys::NP>(double)>;
    // The curved inflow per well, kept so the exact test can evaluate it.
    // The slope change is over the distance to shut-in, floored at 20 bar:
    // a well operating a fraction of a bar below its shut-in would
    // otherwise get a curvature that turns its oil negative a bar away.
    std::vector<double> bhp0s;
    auto curved = [&](Sys& sys, const double kappa, std::vector<Inflow>& inflows) {
        inflows.assign(sys.numWells(), nullptr);
        bhp0s.assign(sys.numWells(), 0.0);
        for (int w = 0; w < sys.numWells(); ++w) {
            const auto& well = sys.wells()[w];
            if (well.shut || well.vfp_table <= 0 || !(well.ipr_b[1] < 0.0)) { continue; }
            const double bhp0 = (well.q_start > 0.0) ? (well.q_start - well.ipr_a[1]) / well.ipr_b[1] : well.bhp_limit;
            bhp0s[w] = bhp0;
            const auto a = well.ipr_a, b = well.ipr_b;
            std::array<double, Sys::NP> c{};
            for (int ph = 0; ph < Sys::NP; ++ph) {
                const double shut_in = (b[ph] < 0.0) ? -a[ph] / b[ph] : bhp0 + 100.0 * unit::barsa;
                c[ph] = kappa * b[ph] / (2.0 * std::max(shut_in - bhp0, 20.0 * unit::barsa));
            }
            inflows[w] = [a, b, c, bhp0](const double bhp) {
                std::array<double, Sys::NP> q{};
                for (int ph = 0; ph < Sys::NP; ++ph) { q[ph] = a[ph] + b[ph] * bhp + c[ph] * (bhp - bhp0) * (bhp - bhp0); }
                return q;
            };
            sys.setWellInflow(w, inflows[w]);
        }
    };
    // The exact decision for the curve: does it miss the tubing at every
    // bhp between the limit and its shut-in? Scanned, since it is a judge.
    auto curve_margin = [&](const Sys& sys, const int w, const Inflow& inflow, const double p_node) {
        const auto& well = sys.wells()[w];
        // The curve's flowing range, found outward from its operating point.
        double lo = bhp0s[w], hi = bhp0s[w];
        while (lo - 0.5 * unit::barsa >= well.bhp_limit && inflow(lo - 0.5 * unit::barsa)[1] > 0.0) { lo -= 0.5 * unit::barsa; }
        while (inflow(hi + 0.5 * unit::barsa)[1] > 0.0 && hi < 1000.0 * unit::barsa) { hi += 0.5 * unit::barsa; }
        double least = std::numeric_limits<double>::max();
        for (int i = 0; i <= 400; ++i) {
            const double bhp = lo + (hi - lo) * i / 400.0;
            auto q = inflow(bhp);
            for (auto& v : q) { v = std::max(v, 0.0); }
            if (!(q[1] > 0.0)) { continue; }
            least = std::min(least, sys.tableBhp(well.vfp_table, p_node, q, well.alq) - well.vfp_dp - bhp);
        }
        return least;
    };
    auto curve_dead = [&](const Sys& sys, const int w, const Inflow& inflow, const double p_node) { return curve_margin(sys, w, inflow, p_node) > 0.0; };
    struct Tally { int n = 0, converged = 0, verified = 0, over_shut = 0, under_shut = 0, wells = 0, shut = 0, pinned_dead = 0,
                   tangent_lifts_curve_not = 0, curve_lifts_tangent_not = 0, margins = 0; long solves = 0, evals = 0;
                   double ms = 0.0, margin_err = 0.0, margin_err_max = 0.0; };
    std::map<std::string, Tally> tally;
    std::vector<std::filesystem::path> files;
    for (const auto& e : std::filesystem::directory_iterator(dir)) { if (e.path().extension() == ".txt") { files.push_back(e.path()); } }
    std::sort(files.begin(), files.end());
    using clock = std::chrono::steady_clock;
    int shown = 0;
    for (const auto& f : files) {
        std::ifstream in(f);
        std::string head; std::getline(in, head);
        if (head != "production") { continue; }
        auto [dumped, guess] = NetworkSolve::readProduction<double>(in, props, units);
        dumped.setAnalyticJacobian(true);
        // The wells the linear solve shuts: given a zero starting rate in the
        // "start shut" runs, so their tangent begins at the bhp limit, as a
        // well the well model has had shut for a while is linearised.
        std::vector<char> shut_by_line(dumped.numWells(), 0);
        {
            auto sys = dumped;
            const auto r = NetworkSolve::solveReduced(sys, guess, params, true);
            if (r.converged) { for (int w = 0; w < sys.numWells(); ++w) { shut_by_line[w] = sys.controlLetter(w) == 'S'; } }
        }
        for (const double kappa : {-0.3, 0.3}) {
            for (const bool rule : {false, true}) {
                for (const bool start_shut : {false, true}) {
                    const bool ext = false;
                    const std::string key = fmt::format("kappa {:+.1f}, tangent {}, start {}", kappa, rule ? "at touch" : "as was ", start_shut ? "shut     " : "as dumped");
                    auto& t = tally[key];
                    auto sys = dumped;
                    if (start_shut) {
                        Sys rebuilt(props, units);
                        for (int nd = 0; nd <= dumped.numNodes(); ++nd) { rebuilt.addNode(dumped.nodes()[nd], dumped.branchAlq(nd)); }
                        rebuilt.setTerminalPressure(dumped.terminalPressure());
                        for (int i = 0; i < dumped.numWells(); ++i) { auto wl = dumped.wells()[i]; if (shut_by_line[i]) { wl.q_start = 0.0; } rebuilt.addWell(wl); }
                        rebuilt.setGroupTarget(dumped.groupTarget());
                        rebuilt.setAnalyticJacobian(true);
                        rebuilt.finish();
                        sys = rebuilt;
                    }
                    std::vector<Inflow> inflows;
                    curved(sys, kappa, inflows);
                    sys.setTangentAtTouchingPoint(rule);
                    const auto t0 = clock::now();
                    std::vector<double> p; std::vector<double> q; bool converged = false; int evals = 0;
                    if (ext) {
                        const auto r = NetworkSolve::solveReducedOnExtension(sys, guess, params);
                        converged = r.converged; p = r.last.node_pressure; q = r.last.well_rate; evals = r.evaluations;
                    } else {
                        const auto r = NetworkSolve::solveReduced(sys, guess, params, true);
                        converged = r.converged; p = r.node_pressure; q = r.well_rate; evals = r.evaluations;
                    }
                    t.ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
                    ++t.n; t.evals += evals; t.solves += sys.wellSolves();
                    if (!converged) { continue; }
                    ++t.converged;
                    const std::string set = letters(sys);
                    if (verifyAnswer(sys, p, q, set).ok) { ++t.verified; }
                    for (int w = 0; w < sys.numWells(); ++w) {
                        if (!inflows[w]) { continue; }
                        const auto& well = sys.wells()[w];
                        const double pn = well.node == 0 ? sys.terminalPressure() : p[well.node];
                        const bool dead = curve_dead(sys, w, inflows[w], pn);
                        if (well.pinned) {
                            // Held as a source by the group; the network does not decide it.
                            if (dead) { ++t.pinned_dead; }
                            continue;
                        }
                        ++t.wells;
                        const bool shut = set[w] == 'S' || !(q[w] > 0.0);
                        if (shut) {
                            ++t.shut;
                            // Where the tangent of a well that cannot lift was
                            // taken: its margin to the tubing against the curve's.
                            // Fractions iterated to the touching bhp, as the rule does.
                            const auto& tbl = props.getTable(well.vfp_table);
                            const double A = detail::getFlo(tbl, well.ipr_a[0], well.ipr_a[1], well.ipr_a[2]);
                            const double B = detail::getFlo(tbl, well.ipr_b[0], well.ipr_b[1], well.ipr_b[2]);
                            double bhp_t = 0.5 * (well.bhp_limit + (well.ipr_b[1] < 0.0 ? -well.ipr_a[1] / well.ipr_b[1] : well.bhp_limit));
                            auto tp = sys.touchingPoint(well, pn, sys.ratesAt(well, bhp_t));
                            for (int pass = 0; pass < 3 && tp.valid && B < 0.0; ++pass) {
                                bhp_t = std::max((tp.flo - A) / B, well.bhp_limit);
                                tp = sys.touchingPoint(well, pn, sys.ratesAt(well, bhp_t));
                            }
                            if (tp.valid && dead) {
                                const double cm = curve_margin(sys, w, inflows[w], pn);
                                const double err = std::abs(tp.gap - cm);
                                ++t.margins; t.margin_err += err; t.margin_err_max = std::max(t.margin_err_max, err);
                                static int outliers = 0;
                                if (trace && rule && err > 5.0 * unit::barsa && outliers++ < 3) {
                                    const double bhp_tan = (well.q_start > 0.0 && well.ipr_b[1] < 0.0) ? (well.q_start - well.ipr_a[1]) / well.ipr_b[1] : well.bhp_limit;
                                    BOOST_TEST_MESSAGE(fmt::format("  outlier {} {}: {} tangent margin {:.2f} curve margin {:.2f} bar at node {:.2f}; touching flo {:.0f} bhp_t {:.1f}; tangent slope {:.0f} sm3/d/bar, line slope {:.0f}; q_start {:.0f} bhp0 {:.1f} limit {:.1f}",
                                                                   f.filename().string(), key, well.name, convert::to(tp.gap, bars), convert::to(cm, bars), convert::to(pn, bars),
                                                                   tp.flo * 86400.0, convert::to(bhp_t, bars), -well.ipr_b[1] * 86400.0 * unit::barsa, -dumped.wells()[w].ipr_b[1] * 86400.0 * unit::barsa,
                                                                   well.q_start * 86400.0, convert::to(bhp0s[w], bars), convert::to(well.bhp_limit, bars)));
                                    (void)bhp_tan;
                                }
                            }
                        }
                        if (shut && !dead) {
                            ++t.over_shut;
                            static int traced = 0;
                            if (trace && rule && start_shut && kappa < 0.0 && traced++ < 3) {
                                const auto& dd = sys.committedDead();
                                const double bhp_tan = well.ipr_b[1] < 0.0 ? -well.ipr_a[1] / well.ipr_b[1] : 0.0;
                                BOOST_TEST_MESSAGE(fmt::format("  sticky {} {}: {} control {} rate {:.0f}, committed dead {}, tangent crossing at settled {:.3f} bar: {:.0f} sm3/d, at guess {:.3f}: {:.0f}; tangent shut-in {:.1f} bar, curve margin {:.2f} bar, q_start {:.0f}",
                                                               f.filename().string(), key, well.name, set[w], q[w] * 86400.0, static_cast<std::size_t>(w) < dd.size() ? int(dd[w]) : -1,
                                                               convert::to(pn, bars), sys.thpPotential(well, pn) * 86400.0, convert::to(guess[well.node], bars), sys.thpPotential(well, guess[well.node]) * 86400.0,
                                                               convert::to(bhp_tan, bars), convert::to(curve_margin(sys, w, inflows[w], pn), bars), well.q_start * 86400.0));
                            }
                        }
                        if (!shut && dead) { ++t.under_shut; }
                        // The IPR question itself: does the final tangent agree
                        // with the curve about lifting at the settled pressure?
                        const double cap = sys.thpPotential(well, pn);
                        const bool tangent_lifts = cap > 0.0;
                        if (tangent_lifts && dead) { ++t.tangent_lifts_curve_not; }
                        if (!tangent_lifts && !dead) {
                            ++t.curve_lifts_tangent_not;
                            if (trace && shown++ < 4) {
                                BOOST_TEST_MESSAGE(fmt::format("  {} {}: {} ({}, {:.0f} sm3/d) tangent says cannot lift at {:.2f} bar, the curve can",
                                                               f.filename().string(), key, well.name, set[w], q[w] * 86400.0, convert::to(pn, bars)));
                            }
                        }
                    }
                }
            }
        }
    }
    for (const auto& [key, t] : tally) {
        BOOST_TEST_MESSAGE(fmt::format("  {}: of {}: converged {}, verified {}; network-decided tubing wells {}, shut {}; shut but the curve lifts {}, open but the curve cannot lift {}; tangent lifts but curve not {}, curve lifts but tangent not {}; pinned wells the curve cannot lift {}; shut wells' margin to the tubing, tangent vs curve: mean {:.3f} max {:.3f} bar over {}; evals/solve {:.1f}, well solves/solve {:.0f}, ms/solve {:.2f}",
                                       key, t.n, t.converged, t.verified, t.wells, t.shut, t.over_shut, t.under_shut, t.tangent_lifts_curve_not, t.curve_lifts_tangent_not, t.pinned_dead,
                                       t.margins ? convert::to(t.margin_err / t.margins, bars) : 0.0, convert::to(t.margin_err_max, bars), t.margins,
                                       t.n ? double(t.evals) / t.n : 0.0, t.n ? double(t.solves) / t.n : 0.0, t.n ? t.ms / t.n : 0.0));
    }
}

// What --network-owns-group-control writes back, judged on dumped systems.
// Every held well must find the group whose own target binds above it, and the
// GRUP targets written to the wells -- each well's share on that group's mode --
// must put the group on its target once the efficiencies are applied up the
// tree. That is the promise the write-back makes to the well model: honour
// these targets and the tree lands where the solve put it.
BOOST_AUTO_TEST_CASE(held_targets_on_the_dumps)
{
    using Sys = DeckTrees::Sys;
    const char* dir = std::getenv("OPM_NETWORK_DUMP_PROD");
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    if (dir == nullptr || inc == nullptr || !std::filesystem::is_directory(dir)) {
        BOOST_TEST_MESSAGE("OPM_NETWORK_DUMP_PROD / OPM_VFP_INCLUDE not set, nothing to check");
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
    const NetworkSolve::Parameters<double> params{1e-2, 80};
    const double r_tol = 0.005;
    auto letters = [](const Sys& sys) {
        std::string s;
        for (int w = 0; w < sys.numWells(); ++w) { s += sys.controlLetter(w); }
        return s;
    };
    struct Tally { int systems = 0, verified = 0, held = 0, unbound = 0, groups_checked = 0, groups_off = 0;
                   std::map<char, int> written; double worst = 0.0; std::string worst_at; };
    std::map<std::string, Tally> tally;
    long wells_total = 0, wells_in_tree = 0;
    int shown = 0;

    for (const auto& e : std::filesystem::directory_iterator(dir)) {
        if (e.path().extension() != ".txt") { continue; }
        std::ifstream in(e.path());
        std::string head; std::getline(in, head);
        if (head != "production") { continue; }
        auto [dumped, guess] = NetworkSolve::readProduction<double>(in, props, units);
        dumped.setAnalyticJacobian(true);
        for (const auto& w : dumped.wells()) { ++wells_total; wells_in_tree += (w.group >= 0); }
        const std::string fname = e.path().filename().string();

        for (const std::string route : {"reduced", "full active set + dead rule"}) {
            auto sys = dumped;
            NetworkSolve::Result<double> r;
            if (route == "reduced") {
                const auto rr = NetworkSolve::solveReduced(sys, guess, params, true);
                r.converged = rr.converged; r.node_pressure = rr.node_pressure; r.well_rate = rr.well_rate;
            } else {
                sys.setDeadWhenCannotLift(true); sys.resetDead();
                r = NetworkSolve::solve(sys, guess, params, NetworkSolve::FullStep{});
            }
            auto& t = tally[route];
            ++t.systems;
            if (!r.converged || !verifyAnswer(sys, r.node_pressure, r.well_rate, letters(sys)).ok) { continue; }
            ++t.verified;
            if (!sys.usesGroupTree()) { continue; }

            const auto held = sys.heldTargets(r.well_rate);
            const auto& groups = sys.groups();
            std::set<int> binding;
            for (int w = 0; w < sys.numWells(); ++w) {
                const char c = sys.controlLetter(w);
                ++t.written[c];
                if (held[w].control != Sys::Control::Tree) { continue; }
                ++t.held;
                if (held[w].group < 0) { ++t.unbound; continue; }
                binding.insert(held[w].group);
            }
            // Each binding group's rate on its own mode, rebuilt from what the
            // write-back hands the wells: a held well's written target, any
            // other well's actual rate, each scaled by the efficiencies above it.
            for (const int g : binding) {
                const auto c = Sys::modeWeights(groups[g].mode, groups[g].resv_coeff);
                double on = 0.0;
                for (int w = 0; w < sys.numWells(); ++w) {
                    const auto& well = sys.wells()[w];
                    double eff = well.efficiency;
                    bool under = false;
                    for (int a = well.group; a >= 0; a = groups[a].parent) {
                        if (a == g) { under = true; break; }
                        eff *= groups[a].efficiency;
                    }
                    if (!under) { continue; }
                    double rate = 0.0;
                    if (held[w].control == Sys::Control::Tree && held[w].group == g) {
                        rate = held[w].value;
                    } else if (well.ipr_b[1] < 0.0) {
                        const double bhp = (r.well_rate[w] - well.ipr_a[1]) / well.ipr_b[1];
                        for (int ph = 0; ph < Sys::NP; ++ph) {
                            rate += c[ph] * std::max(well.ipr_a[ph] + well.ipr_b[ph] * bhp, 0.0);
                        }
                    }
                    on += eff * rate;
                }
                ++t.groups_checked;
                // As the judge: relative, but never tighter than 1 sm3/d, so a
                // zero target is not violated by rounding.
                const double miss = std::abs(on - groups[g].target);
                const double off = miss / std::max(groups[g].target, 1.0 / 86400.0);
                if (off > t.worst) { t.worst = off; t.worst_at = fname + " " + groups[g].name; }
                if (miss > std::max(r_tol * groups[g].target, 1.0 / 86400.0)) {
                    ++t.groups_off;
                    if (shown++ < 6) {
                        BOOST_TEST_MESSAGE(fmt::format("  {} {}: {} on its mode from the written targets is {:.4f} of its target",
                                                       route, fname, groups[g].name, on / groups[g].target));
                    }
                }
            }
        }
    }
    BOOST_TEST_MESSAGE(fmt::format("wells in the tree: {} of {} across the dumps", wells_in_tree, wells_total));
    for (const auto& [route, t] : tally) {
        std::string w;
        for (const auto& [c, n] : t.written) { w += fmt::format(" {}x{}", c, n); }
        BOOST_TEST_MESSAGE(fmt::format("  {:28} {} systems, {} verified; held {} ({} with no binding group); "
                                       "binding groups checked {}, off target {} (worst {:.2e} at {}); controls{}",
                                       route, t.systems, t.verified, t.held, t.unbound, t.groups_checked,
                                       t.groups_off, t.worst, t.worst_at.empty() ? "-" : t.worst_at, w));
        BOOST_CHECK_EQUAL(t.unbound, 0);
        BOOST_CHECK_EQUAL(t.groups_off, 0);
    }
}

// The same systems solved twice, the IPR as given: once from the guessed
// pressures, once from Stein's allocation and the pressures it implies. Stein's
// balancer needs a GuideRate, hence a Schedule: OPM_NETWORK_DECK names the deck
// the dumps came from. Also asks whether each answer is the allocation Stein's
// balancer gives at the answer's own pressures -- the same tree, GCONPROD item 8
// included.
BOOST_AUTO_TEST_CASE(stein_start_on_the_dumps)
{
    using Sys = DeckTrees::Sys;
    const char* dir = std::getenv("OPM_NETWORK_DUMP_PROD");
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    const char* deck_path = std::getenv("OPM_NETWORK_DECK");
    if (dir == nullptr || inc == nullptr || deck_path == nullptr || !std::filesystem::is_directory(dir)) {
        BOOST_TEST_MESSAGE("OPM_NETWORK_DUMP_PROD / OPM_VFP_INCLUDE / OPM_NETWORK_DECK not set, nothing to run");
        return;
    }
    std::deque<VFPProdTable> tables;
    VFPProdProperties<double> props;
    const UnitSystem units{};
    for (const char* name : {"well_vfp.ecl", "flowl_b_vfp.ecl", "flowl_c_vfp.ecl"}) {
        const auto path = std::filesystem::path(inc) / name;
        if (!std::filesystem::exists(path)) { continue; }
        const auto d = Parser{}.parseFile(path.string());
        for (const auto& kw : d.getKeywordList("VFPPROD")) {
            tables.emplace_back(*kw, /*gaslift_opt_active=*/true, units);
            props.addTable(tables.back());
        }
    }
    const auto deck = Parser{}.parseFile(deck_path);
    const EclipseState es(deck);
    const Schedule schedule(deck, es);

    const NetworkSolve::Parameters<double> params{1e-2, 80};
    auto letters = [](const Sys& sys) {
        std::string s;
        for (int w = 0; w < sys.numWells(); ++w) { s += sys.controlLetter(w); }
        return s;
    };
    struct Tally { int n = 0, converged = 0, verified = 0, stein_agrees = 0; long iterations = 0, set_changes = 0; };
    std::map<std::string, Tally> tally;
    int stein_failed = 0, systems = 0;
    std::map<std::string, std::string> flips;      // system -> what the start changed

    auto trace_changes = [](const std::string& trace) {
        long k = 0; std::string prev;
        std::istringstream in(trace);
        for (std::string t; in >> t; ) { if (!prev.empty() && t != prev) { ++k; } prev = t; }
        return k;
    };

    for (const auto& e : std::filesystem::directory_iterator(dir)) {
        if (e.path().extension() != ".txt") { continue; }
        std::ifstream in(e.path());
        std::string head; std::getline(in, head);
        if (head != "production") { continue; }
        auto [dumped, guess] = NetworkSolve::readProduction<double>(in, props, units);
        dumped.setAnalyticJacobian(true);
        if (!dumped.usesGroupTree()) { continue; }
        ++systems;
        const std::string fname = e.path().filename().string();

        // The start, computed once on a copy and handed to both routes.
        auto probe = dumped;
        GuideRate gr_start{schedule};
        const auto start = NetworkSolve::steinStart(probe, guess, gr_start, 0);
        if (!start.ok) { ++stein_failed; }

        for (const std::string route : {"reduced", "full active set + dead rule"}) {
            for (const bool from_stein : {false, true}) {
                if (from_stein && !start.ok) { continue; }
                auto sys = dumped;
                const auto& p0 = from_stein ? start.node_pressure : guess;
                bool conv = false; std::vector<double> p, q; long it = 0, changes = 0;
                if (route == "reduced") {
                    const auto r = NetworkSolve::solveReduced(sys, p0, params, true);
                    conv = r.converged; p = r.node_pressure; q = r.well_rate; it = r.iterations; changes = r.set_changes;
                } else {
                    sys.setDeadWhenCannotLift(true); sys.resetDead();
                    if (from_stein) { sys.setStartAllocation(start.well_rate); }
                    const auto r = NetworkSolve::solve(sys, p0, params, NetworkSolve::FullStep{});
                    conv = r.converged; p = r.node_pressure; q = r.well_rate; it = r.iterations;
                    changes = trace_changes(r.control_trace);
                }
                const std::string key = route + (from_stein ? "  from Stein" : "  from guess");
                auto& t = tally[key];
                ++t.n;
                if (!conv) {
                    if (from_stein) { flips[fname] += fmt::format(" [{}: guess {} -> Stein fail]", route, "?"); }
                    continue;
                }
                ++t.converged; t.iterations += it; t.set_changes += changes;
                const bool ok = verifyAnswer(sys, p, q, letters(sys)).ok;
                t.verified += ok;
                // Stein's allocation at this answer's pressures.
                auto oracle_sys = dumped;
                GuideRate gr_oracle{schedule};
                const auto x = oracle_sys.start(p);
                const auto stein = NetworkSolve::steinAllocation(oracle_sys, x, gr_oracle, 0);
                if (!stein.empty()) {
                    double worst = 0.0;
                    for (int w = 0; w < sys.numWells(); ++w) {
                        worst = std::max(worst, std::abs(q[w] - stein[w]) / std::max(std::abs(stein[w]), 1.0 / 86400.0));
                    }
                    t.stein_agrees += (worst < 0.005);
                    static int shown = 0;
                    if (worst >= 0.005 && shown < 4 && route == "reduced" && !from_stein) {
                        ++shown;
                        std::string line;
                        for (int w = 0; w < sys.numWells(); ++w) {
                            line += fmt::format(" {}={:.0f}/{:.0f}({})", sys.wells()[w].name, q[w] * 86400.0,
                                                stein[w] * 86400.0, sys.controlLetter(w));
                        }
                        std::string binds = sys.treeSignature();
                        BOOST_TEST_MESSAGE(fmt::format("  {} ours/Stein oil sm3/d:{}  set {}", fname, line, binds));
                    }
                }
            }
        }
    }
    BOOST_TEST_MESSAGE(fmt::format("{} systems with a tree; Stein's balancer rejected {} of them", systems, stein_failed));
    for (const auto& [key, t] : tally) {
        BOOST_TEST_MESSAGE(fmt::format("  {:44} of {}: converged {}, verified {}, same as Stein's allocation {}, "
                                       "iterations/solve {:.2f}, set changes/solve {:.2f}",
                                       key, t.n, t.converged, t.verified, t.stein_agrees,
                                       t.converged ? double(t.iterations) / t.converged : 0.0,
                                       t.converged ? double(t.set_changes) / t.converged : 0.0));
    }
}

BOOST_AUTO_TEST_SUITE_END()
