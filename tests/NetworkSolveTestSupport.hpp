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
#ifndef OPM_NETWORK_SOLVE_TEST_SUPPORT_HEADER_INCLUDED
#define OPM_NETWORK_SOLVE_TEST_SUPPORT_HEADER_INCLUDED

// What the network bench cases are built from: the deck and reference
// fixtures, the standalone Newton and its globalisations, the generated
// trees and decks, and the judge. Shared by the four case files, so each
// of those is only its own subject.

/*!
 * \file
 *
 * \brief A standalone bench for the injection-network solve.
 *
 * The GNETINJE_GAS-01 network with the wells replaced by their inflow
 * performance plus the control logic, so solution methods can be compared
 * without a reservoir or a well solve. The VFPINJ tables are the deck's.
 *
 * Topology (GNETINJE_GAS-01, table 9999 = no table, pressure passes through):
 *
 *     PLAT-A  340 bar terminal
 *       |  VFPINJ 3          total gas
 *      M5S ------------------------------- G1 (9999)  G-3H, G-4H
 *       |  VFPINJ 2          F-wells' gas
 *      M5N ------------------------------- F1 (9999)  F-1H, F-2H
 *
 * so there are two unknowns, p(M5S) and p(M5N), and
 *
 *     G(p)_M5S = vfp3.bhp(thp = 340 bar, q = sum of all four wells)
 *     G(p)_M5N = vfp2.bhp(thp = p_M5S,   q = q(F-1H) + q(F-2H))
 *
 * with each well's rate taken from its own VFPINJ 1 against a linear IPR and
 * then put through the control logic (THP / BHP limit / rate limit / group
 * target). That last part is what gives the response its plateau, and it is
 * the part a pure dq/dbhp proxy misses.
 *
 * The bench is calibrated to the reference solution at day 31 and reproduces
 * it, so the numbers it reports are about the methods and not about the model.
 * What it is mainly for is globalisation: the Newton direction is the same in
 * FullStep, CappedStep, LineSearch and TrustRegion, and globalisation_basin
 * measures how much of the starting-pressure space each one recovers from.
 */

#include <config.h>

#define BOOST_TEST_MODULE NetworkSolveBench

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/input/eclipse/Deck/UDAValue.hpp>
#include <opm/io/eclipse/ESmry.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/simulators/wells/ProdGroupTreeNode.hpp>
#include <opm/simulators/wells/ProdGroupTreeBalancer.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/Network/Branch.hpp>
#include <opm/input/eclipse/Schedule/Network/Node.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>
#include <opm/input/eclipse/Schedule/VFPProdTable.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/simulators/wells/network/NetworkNodePressureUpdater.hpp>
#include <opm/simulators/wells/network/NetworkInjectionSystem.hpp>
#include <opm/simulators/wells/network/NetworkJudge.hpp>
#include <opm/simulators/wells/network/NetworkProductionSystem.hpp>
#include <opm/simulators/wells/network/NetworkLegacySolve.hpp>
#include <opm/simulators/wells/network/NetworkReducedSolve.hpp>
#include <opm/simulators/wells/network/NetworkTreeSolve.hpp>
#include <opm/simulators/wells/network/NetworkTubingExtension.hpp>

#include <algorithm>
#include <array>
#include <deque>
#include <set>
#include <filesystem>
#include <sstream>
#include <fmt/format.h>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <tuple>
#include <map>
#include <memory>
#include <iomanip>
#include <cmath>
#include <random>
#include <chrono>
#include <functional>
#include <string>
#include <vector>

using namespace Opm;
using namespace Opm::unit;

namespace {

/// The settings every bench case uses unless it is varying them on purpose.
constexpr NetworkSolve::Parameters<double> kParams{1e-2, 50};

// VFPINJ 1 (wells), from opm-tests/network/include/vfp_gi_wells.inc.
const std::string vfp_well = R"(
VFPINJ
      1        2011         GAS /
-- gas rates Sm3/d
    5000   42642   92830  168113  243396  368868  494340  619811
  745283  870755  996226 1121698 1247170 1498113 1749057 2000000 /

-- Tubing head pressure [bar]
  50.00 100.00 150.00 200.00 250.00 300.00 350.00 400.00
 450.00 500.00 /

  1   71.916  71.436  69.906  65.520  58.044  32.523   0.000
   0.000   0.000   0.000   0.000   0.000   0.000   0.000
   0.000   0.000 /

  2  148.936 148.783 148.181 146.416 143.591 136.308 125.221
 109.013  84.309  33.054   0.000   0.000   0.000   0.000
   0.000   0.000 /

  3  223.875 223.824 223.518 222.519 220.866 216.643 210.493
 202.241 191.582 178.067 160.870 138.369 106.524   0.000
   0.000   0.000 /

  4  292.521 292.511 292.307 291.583 290.369 287.238 282.709
 276.721 269.214 260.055 249.100 236.105 220.723 180.046
 111.512   0.000 /

  5  356.485 356.465 356.302 355.710 354.690 352.069 348.274
 343.287 337.085 329.619 320.806 310.596 298.855 270.132
 232.423 180.015 /

  6  417.573 417.471 417.328 416.798 415.900 413.565 410.188
 405.762 400.274 393.685 385.994 377.130 367.053 342.950
 312.850 275.232 /

  7  476.662 476.560 476.427 475.948 475.111 472.959 469.858
 465.789 460.760 454.742 447.714 439.687 430.578 409.026
 382.638 350.712 /

  8  534.424 534.322 534.200 533.741 532.966 530.946 528.029
 524.214 519.502 513.882 507.323 499.836 491.401 471.501
 447.388 418.644 /

  9  591.218 591.126 591.004 590.565 589.821 587.903 585.129
 581.508 577.030 571.695 565.494 558.415 550.438 531.742
 509.190 482.537 /

 10  647.287 647.185 647.063 646.645 645.931 644.085 641.422
 637.954 633.660 628.560 622.624 615.851 608.242 590.443
 569.053 543.900 /
)";

// VFPINJ 2 (M5S -> M5N) and 3 (PLAT-A -> M5S), same source directory.
const std::string vfp_m5n = R"(
VFPINJ
      2        288.0         GAS /
-- gas rates Sm3/d
    5000   35227   85606  135985  186364  236742  287121  337500
  387879  488636  589394  690151  790909  891667  992424 1193939
 1395454 1596969 1798485 2000000 /

-- Tubing head pressure [bar]
  50.00 100.00 150.00 200.00 250.00 300.00 350.00 400.00
 450.00 500.00 /

  1   56.550  55.729  52.618  46.743  36.548  12.637   0.000
   0.000   0.000   0.000   0.000   0.000   0.000   0.000
   0.000   0.000   0.000   0.000   0.000   0.000 /

  2  112.720 112.299 110.960 108.660 105.290 100.711  94.706
  86.919  76.702  40.155   0.000   0.000   0.000   0.000
   0.000   0.000   0.000   0.000   0.000   0.000 /

  3  168.848 168.470 167.584 166.126 164.042 161.299 157.864
 153.685 148.684 135.832 117.872  91.120  34.507   0.000
   0.000   0.000   0.000   0.000   0.000   0.000 /

  4  224.533 224.166 223.453 222.297 220.688 218.593 216.001
 212.890 209.251 200.254 188.731 174.172 155.694 131.383
  95.743   0.000   0.000   0.000   0.000   0.000 /

  5  279.786 279.451 278.814 277.820 276.438 274.645 272.442
 269.828 266.772 259.330 250.021 238.659 225.008 208.689
 189.033 133.046   0.000   0.000   0.000   0.000 /

  6  334.758 334.434 333.861 332.954 331.712 330.103 328.126
 325.783 323.061 316.473 308.287 298.437 286.816 273.252
 257.527 218.085 161.644   0.000   0.000   0.000 /

  7  389.535 389.233 388.693 387.861 386.706 385.215 383.390
 381.230 378.724 372.666 365.192 356.250 345.774 333.678
 319.854 286.363 243.130 184.314  71.043   0.000 /

  8  444.183 443.902 443.384 442.596 441.505 440.101 438.384
 436.353 434.010 428.340 421.352 413.025 403.305 392.149
 379.480 349.273 311.592 264.320 202.058  94.587 /

  9  498.745 498.464 497.978 497.222 496.174 494.846 493.215
 491.282 489.046 483.657 477.037 469.164 459.994 449.497
 437.638 409.580 375.150 333.246 281.740 215.288 /

 10  553.220 552.961 552.486 551.762 550.758 549.472 547.906
 546.049 543.910 538.748 532.408 524.881 516.133 506.143
 494.878 468.364 436.137 397.516 351.346 295.402 /
)";

const std::string vfp_m5s = R"(
VFPINJ
      3        285.0         GAS /
-- gas rates Sm3/d
    5000   35227   85606  135985  186364  236742  287121  337500
  387879  488636  589394  690151  790909  891667  992424 1193939
 1395454 1596969 1798485 2000000 /

-- Tubing head pressure [bar]
  50.00 100.00 150.00 200.00 250.00 300.00 350.00 400.00
 450.00 500.00 /

  1   68.834  67.861  64.174  57.211  45.128  16.789   0.000
   0.000   0.000   0.000   0.000   0.000   0.000   0.000
   0.000   0.000   0.000   0.000   0.000   0.000 /

  2  135.406 134.907 133.320 130.594 126.600 121.173 114.056
 104.827  92.718  49.403   0.000   0.000   0.000   0.000
   0.000   0.000   0.000   0.000   0.000   0.000 /

  3  201.928 201.480 200.430 198.702 196.232 192.981 188.910
 183.957 178.030 162.798 141.512 109.806  42.709   0.000
   0.000   0.000   0.000   0.000   0.000   0.000 /

  4  267.925 267.490 266.645 265.275 263.368 260.885 257.813
 254.126 249.813 239.150 225.493 208.238 186.338 157.525
 115.285   0.000   0.000   0.000   0.000   0.000 /

  5  333.410 333.013 332.258 331.080 329.442 327.317 324.706
 321.608 317.986 309.166 298.133 284.667 268.488 249.147
 225.851 159.496   0.000   0.000   0.000   0.000 /

  6  398.562 398.178 397.499 396.424 394.952 393.045 390.702
 387.925 384.699 376.891 367.189 355.515 341.742 325.666
 307.029 260.283 193.390   0.000   0.000   0.000 /

  7  463.483 463.125 462.485 461.499 460.130 458.363 456.200
 453.640 450.670 443.490 434.632 424.034 411.618 397.282
 380.898 341.205 289.966 220.258  86.011   0.000 /

  8  528.251 527.918 527.304 526.370 525.077 523.413 521.378
 518.971 516.194 509.474 501.192 491.323 479.803 466.581
 451.566 415.765 371.106 315.080 241.288 113.915 /

  9  592.917 592.584 592.008 591.112 589.870 588.296 586.363
 584.072 581.422 575.035 567.189 557.858 546.990 534.549
 520.494 487.240 446.434 396.770 335.726 256.968 /

 10  657.480 657.173 656.610 655.752 654.562 653.038 651.182
 648.981 646.446 640.328 632.814 623.893 613.525 601.685
 588.334 556.910 518.715 472.942 418.222 351.918 /
)";


/// A VFPPROD table, from tests/test_networkpressure.cpp: LIQ rate, WCT / GOR
/// fractions, GRAT alq. Small enough to reason about by hand.
const std::string vfp_prod = R"(
VFPPROD
     3     250.00      LIQ        WCT         GOR         THP        GRAT      METRIC   BHP      /
       20.0       100.0    1000.0     2000.0 /
      10.00      30.00 /
      0.000      0.5      1.0 /
       100.0 /
        0.0 /
  1  1  1  1    12.0   15.0   20.0   30.0 /
  1  2  1  1    13.0   16.0   21.0   31.0 /
  1  3  1  1    14.0   17.0   22.0   32.0 /
  2  1  1  1    32.0   35.0   40.0   50.0 /
  2  2  1  1    33.0   36.0   41.0   51.0 /
  2  3  1  1    34.0   37.0   42.0   52.0 /
)";

// ---------------------------------------------------------------------------
// Model
//
// A network case is data: nodes with a parent and a VFP table, wells hanging
// off nodes, a terminal pressure and an optional group target. Nothing below
// knows about GNETINJE_GAS-01 in particular -- gnetinjeGas() is one instance,
// and NetworkCase::Builder is what a deck reader would fill in.
// ---------------------------------------------------------------------------

constexpr int kNoTable = 9999;   // GNETINJE's "no table": pressure passes through

/// Which phase the network carries. A production network would need VFPPROD
/// instead, and with it a water and a gas fraction per branch and an ALQ -- see
/// the note on Rates below.
enum class Fluid { Gas, Water };

/// The rate triple a VFP lookup takes. For an injection network only one entry
/// is ever non-zero, but carrying the triple is what a production network would
/// need: there the branch rate splits into oil, water and gas, and VFPPROD is
/// looked up on a flow rate plus WFR and GFR fractions (and an ALQ). Those
/// fractions are extra unknowns per branch, with their own mixing equations at
/// the nodes -- which is why the production side is a bigger job than swapping
/// the table type.
struct Rates
{
    double aqua = 0.0;
    double liquid = 0.0;
    double vapour = 0.0;
};

class Reference;

/// A network node. Node 0 is the terminal and carries the fixed pressure.
struct Node
{
    std::string name;
    int parent = -1;        // -1 only for the terminal
    int vfp_table = kNoTable;
};

/// One injector: a linear IPR against its own tubing table, plus its limits.
struct Well
{
    std::string name;
    int node = 0;
    int vfp_table = 1;
    double q_ref = 0.0;      // rate in the reference solution            [sm3/s]
    double bhp_ref = 0.0;    // its own bhp there; the IPR pivots here       [Pa]
    double dq_dbhp = 0.0;    // IPR slope, the stiffness knob        [sm3/s/Pa]
    double bhp_limit = 0.0;
    double rate_limit = 0.0;
    double guide = 0.0;      // share of a group target; defaults to q_ref
    double efficiency = 1.0; // WEFAC as the network sees it
};

class NetworkCase
{
public:
    /// Add a VFPINJ table from deck text. The table number in the text is the
    /// one the nodes and wells refer to.
    void addTable(const std::string& deck_text)
    {
        decks_.push_back(Parser{}.parseString(deck_text));
        addInjTable(decks_.back()["VFPINJ"].front());
    }

    void addInjTable(const DeckKeyword& keyword)
    {
        tables_.emplace_back(keyword, UnitSystem{});
        props_.addTable(tables_.back());

        const auto& t = tables_.back();
        axes_[t.getTableNum()] = Axes{t.getFloAxis().front(), t.getFloAxis().back(),
                                      t.getTHPAxis().front(), t.getTHPAxis().back()};
    }

    void setFluid(const Fluid f) { fluid_ = f; }
    Fluid fluid() const { return fluid_; }

    /// Apply an operating point: this is what fixes each well's bhp_ref, and so
    /// what its IPR is a linearisation about.
    void calibrate(const Reference& reference);

    void addNode(Node n) { nodes_.push_back(std::move(n)); }
    void addWell(Well w) { wells_.push_back(std::move(w)); }
    void setTerminalPressure(const double p) { terminal_pressure_ = p; }
    void setGroupTarget(const double target) { group_target_ = target; }

    /// Resolve everything derived from the reference solution. Call once the
    /// nodes, wells and tables are in.
    void finish()
    {
        for (auto& w : wells_) {
            if (w.guide <= 0.0) {
                w.guide = w.q_ref;
            }
        }
        children_.assign(nodes_.size(), {});
        wells_at_.assign(nodes_.size(), {});
        for (std::size_t n = 1; n < nodes_.size(); ++n) {
            children_[nodes_[n].parent].push_back(static_cast<int>(n));
        }
        for (std::size_t w = 0; w < wells_.size(); ++w) {
            wells_at_[wells_[w].node].push_back(static_cast<int>(w));
        }
        // Nodes whose pressure is not simply their parent's: the real unknowns
        // of the eliminated form.
        solved_.clear();
        for (std::size_t n = 1; n < nodes_.size(); ++n) {
            if (hasTable(nodes_[n])) {
                solved_.push_back(static_cast<int>(n));
            }
        }
    }

    /// IPR slope shared by all wells [sm3/d per bar]: the knob the bench exists for.
    void setStiffness(const double dq_dbhp_sm3_day_per_bar)
    {
        const double si = convert::from(dq_dbhp_sm3_day_per_bar, cubic(meter) / day)
                        / convert::from(1.0, bars);
        for (auto& w : wells_) {
            w.dq_dbhp = si;
        }
    }

    /// The same knob as a fraction of each well's own rate per bar, which is the
    /// only form that transfers between cases: a gas injector taking 5e5 sm3/d
    /// and a water injector taking 700 have nothing comparable to say in
    /// absolute units. 0.12/bar reproduces the 6e4 sm3/d/bar used for the gas
    /// case. Call after the rates are known.
    void setRelativeStiffness(const double fraction_per_bar)
    {
        for (auto& w : wells_) {
            w.dq_dbhp = fraction_per_bar * w.q_ref / convert::from(1.0, bars);
        }
    }

    /// Clamp table lookups to the flow and THP axes, as the simulator's network
    /// pressure computation does. See table_bounds_want_to_be_constraints.
    void setClampToAxes(const bool on) { clamp_to_axes_ = on; }

    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<Well>& wells() const { return wells_; }
    std::vector<Well>& wells() { return wells_; }
    const std::vector<int>& children(const int n) const { return children_[n]; }
    const std::vector<int>& wellsAt(const int n) const { return wells_at_[n]; }
    const std::vector<int>& solvedNodes() const { return solved_; }
    double terminalPressure() const { return terminal_pressure_; }
    double groupTarget() const { return group_target_; }
    bool hasTable(const Node& n) const { return n.vfp_table != kNoTable && axes_.count(n.vfp_table); }

    /// The library system this case describes: the same object the simulator
    /// assembles, differing only in where the wells' inflow performance came
    /// from. Here it is a linearisation about the reference operating point.
    NetworkSolve::InjectionSystem<double> system() const
    {
        NetworkSolve::InjectionSystem<double> s(props_, fluid_ == Fluid::Gas ? Phase::GAS : Phase::WATER);
        s.setTerminalPressure(terminal_pressure_);
        s.setGroupTarget(group_target_);
        s.setClampToAxes(clamp_to_axes_);
        for (const auto& n : nodes_) {
            s.addNode(NetworkSolve::Node{n.name, n.parent, n.vfp_table});
        }
        for (const auto& w : wells_) {
            NetworkSolve::Well<double> sw;
            sw.name = w.name;
            sw.node = w.node;
            sw.vfp_table = w.vfp_table;
            // q = q_ref + dq_dbhp*(bhp - bhp_ref) as q = a + b*bhp.
            sw.ipr_a = w.dq_dbhp * w.bhp_ref - w.q_ref;
            sw.ipr_b = w.dq_dbhp;
            sw.bhp_limit = w.bhp_limit;
            sw.rate_limit = w.rate_limit;
            sw.guide = w.guide > 0.0 ? w.guide : w.q_ref;
            sw.efficiency = w.efficiency;
            sw.in_group = group_target_ > 0.0;
            s.addWell(sw);
        }
        s.finish();
        return s;
    }

    /// Rebuild a system written by the simulator, against this case's tables.
    std::pair<NetworkSolve::InjectionSystem<double>, std::vector<double>>
    systemFromDump(std::istream& is) const
    {
        return NetworkSolve::read<double>(is, props_);
    }

    /// The network's scalar rate as the triple a VFP lookup takes.
    Rates asRates(const double q) const
    {
        Rates r;
        (fluid_ == Fluid::Gas ? r.vapour : r.aqua) = q;
        return r;
    }

    /// Downstream pressure of a branch, or a well's bhp: the same table lookup.
    double tableBhp(const int table, const double thp, const double q) const
    {
        const auto r = asRates(clamp_to_axes_
            ? std::clamp(q, axes_.at(table).flo_min, axes_.at(table).flo_max) : q);
        const double p = clamp_to_axes_
            ? std::clamp(thp, axes_.at(table).thp_min, axes_.at(table).thp_max) : thp;
        return props_.bhp(table, r.aqua, r.liquid, r.vapour, p);
    }

    /// Largest rate the table describes. Past it the cells are zero-filled and
    /// the interpolation runs away, so this is the edge of the feasible set.
    double maxFlow(const int table) const { return axes_.at(table).flo_max; }

    static double ipr(const Well& w, const double bhp)
    {
        return w.q_ref + w.dq_dbhp * (bhp - w.bhp_ref);
    }

    /// The rate this well takes at a given node pressure, with its own limits
    /// applied -- the eliminated form's inner solve.
    double wellRate(const Well& w, const double p_node) const
    {
        const auto f = [&](const double q) { return ipr(w, tableBhp(w.vfp_table, p_node, q)) - q; };
        double lo = convert::from(5000.0, cubic(meter) / day);
        double hi = w.rate_limit;
        while (hi > lo && tableBhp(w.vfp_table, p_node, hi) <= convert::from(1.0, atm)) {
            hi *= 0.9;
        }
        if (f(lo) <= 0.0) {
            return 0.0;                       // IPR cannot deliver even the axis minimum
        }
        double q = hi;
        if (f(hi) < 0.0) {
            for (int it = 0; it < 60; ++it) {
                q = 0.5 * (lo + hi);
                (f(q) > 0.0 ? lo : hi) = q;
            }
        }
        if (tableBhp(w.vfp_table, p_node, q) > w.bhp_limit) {
            q = std::max(ipr(w, w.bhp_limit), 0.0);   // BHP-limited
        }
        return std::clamp(q, 0.0, w.rate_limit);
    }

    /// Per-well rates at the given node pressures, group target applied.
    std::vector<double> rates(const std::vector<double>& node_pressure) const
    {
        std::vector<double> q(wells_.size());
        for (std::size_t i = 0; i < wells_.size(); ++i) {
            q[i] = wellRate(wells_[i], node_pressure[wells_[i].node]);
        }
        if (group_target_ > 0.0) {
            const double sum = std::accumulate(q.begin(), q.end(), 0.0);
            if (sum > group_target_) {
                // GRUP control: share the target out by guide rate.
                double guides = 0.0;
                for (const auto& w : wells_) {
                    guides += w.guide;
                }
                for (std::size_t i = 0; i < q.size(); ++i) {
                    q[i] = std::min(q[i], group_target_ * wells_[i].guide / guides);
                }
            }
        }
        return q;
    }

    /// Pressure at every node, given the pressures applied to the solved ones.
    std::vector<double> nodePressures(const std::vector<double>& applied) const
    {
        std::vector<double> p(nodes_.size(), terminal_pressure_);
        for (std::size_t n = 1; n < nodes_.size(); ++n) {
            const auto it = std::find(solved_.begin(), solved_.end(), static_cast<int>(n));
            p[n] = (it != solved_.end()) ? applied[it - solved_.begin()] : p[nodes_[n].parent];
        }
        return p;
    }

    /// Rate through each node's parent branch, from the well rates upwards.
    std::vector<double> branchFlows(const std::vector<double>& well_rate) const
    {
        std::vector<double> q(nodes_.size(), 0.0);
        for (std::size_t n = nodes_.size(); n-- > 1;) {
            for (const int w : wells_at_[n]) {
                q[n] += well_rate[w];
            }
            for (const int c : children_[n]) {
                q[n] += q[c];
            }
        }
        return q;
    }

private:
    struct Axes { double flo_min, flo_max, thp_min, thp_max; };

    std::map<int, Axes> axes_;
    std::deque<Deck> decks_;
    std::deque<VFPInjTable> tables_;
    VFPInjProperties<double> props_;

    std::vector<Node> nodes_;
    std::vector<Well> wells_;
    std::vector<std::vector<int>> children_;
    std::vector<std::vector<int>> wells_at_;
    std::vector<int> solved_;

    Fluid fluid_ = Fluid::Gas;
    double terminal_pressure_ = 0.0;
    double group_target_ = 0.0;
    bool clamp_to_axes_ = false;
};

/// The operating point a case is calibrated against: for each well, the rate it
/// takes and its own bottom-hole pressure. Set it by hand -- which is what the
/// unit tests do, and what you would do to pose a difficult case -- or read it
/// from a reference run's summary.
///
/// It has to be the well's bhp and not the pressure at its node. The two agree
/// only while the well is on THP control; under group control the well sits well
/// below its node (162 bar below, at day 91 of GNETINJE_GAS-01), and calibrating
/// against the node pressure there produces a well that does not exist.
class Reference
{
public:
    /// Rate in sm3/d, bottom-hole pressure in bar.
    void set(const std::string& well, const double rate_sm3_day, const double bhp_bar)
    {
        point_[well] = {convert::from(rate_sm3_day, cubic(meter) / day),
                        convert::from(bhp_bar, bars)};
    }

    bool has(const std::string& well) const { return point_.count(well) > 0; }
    double rate(const std::string& well) const { return point_.at(well).first; }
    double bhp(const std::string& well) const { return point_.at(well).second; }

    /// Read the point out of a summary case at the given time, taking the last
    /// report at or before it. `prefix` is the case path without .SMSPEC.
    /// Throws if a vector is missing, so a caller that wants to skip when the
    /// reference run is not available should check the file exists first.
    static Reference fromSummary(const std::string& prefix,
                                 const NetworkCase& c,
                                 const double time_days);

private:
    std::map<std::string, std::pair<double, double>> point_;
};

/// Topology, tables and limits from a deck; nothing calibrated yet. Reads the
/// keywords straight off the parsed Deck rather than building a Schedule, so it
/// stays usable from a unit test. WELSPECS and GRUPTREE are accumulated over
/// every occurrence; GNETINJE and WCONINJE are taken from their first, which is
/// the state at the start of the run rather than at any later DATES.
NetworkCase fromDeck(const std::string& deck_path, const Fluid fluid)
{
    const auto deck = Parser{}.parseFile(deck_path);
    NetworkCase c;
    c.setFluid(fluid);

    // GNETINJE spells the phase WAT where WCONINJE spells it WATER.
    const std::string well_phase = (fluid == Fluid::Gas) ? "GAS" : "WATER";
    const std::string node_phase = (fluid == Fluid::Gas) ? "GAS" : "WAT";

    // Which group each well belongs to, and each group's parent. These can be
    // spread over several keywords, so take them all.
    std::map<std::string, std::string> well_group;
    for (const auto& keyword : deck["WELSPECS"]) {
        for (const auto& record : keyword) {
            well_group[record.getItem("WELL").getTrimmedString(0)] =
                record.getItem("GROUP").getTrimmedString(0);
        }
    }
    std::map<std::string, std::string> parent;
    for (const auto& keyword : deck["GRUPTREE"]) {
        for (const auto& record : keyword) {
            parent[record.getItem("CHILD_GROUP").getTrimmedString(0)] =
                record.getItem("PARENT_GROUP").getTrimmedString(0);
        }
    }

    // The network itself: which groups are nodes, their table, and the terminal.
    std::map<std::string, int> node_table;
    std::string terminal;
    for (const auto& record : deck["GNETINJE"].front()) {
        if (record.getItem("PHASE").getTrimmedString(0) != node_phase) {
            continue;
        }
        const auto name = record.getItem("GROUP").getTrimmedString(0);
        const auto& pressure = record.getItem("PRESSURE");
        if (pressure.hasValue(0) && !pressure.defaultApplied(0)) {
            terminal = name;
            c.setTerminalPressure(pressure.getSIDouble(0));
        }
        const auto& table = record.getItem("VFP_TABLE");
        node_table[name] = (table.hasValue(0) && !table.defaultApplied(0))
            ? table.get<int>(0) : kNoTable;
    }
    if (terminal.empty()) {
        throw std::runtime_error("no terminal node in GNETINJE for " + node_phase);
    }

    // Nodes, terminal first, then each node after its parent.
    std::vector<std::string> order{terminal};
    for (bool grew = true; grew;) {
        grew = false;
        for (const auto& [name, table] : node_table) {
            const bool placed = std::find(order.begin(), order.end(), name) != order.end();
            const auto p = parent.find(name);
            if (placed || p == parent.end()) {
                continue;
            }
            const auto at = std::find(order.begin(), order.end(), p->second);
            if (at != order.end()) {
                order.push_back(name);
                grew = true;
            }
        }
    }
    auto index = [&order](const std::string& name) {
        return static_cast<int>(std::find(order.begin(), order.end(), name) - order.begin());
    };
    for (std::size_t i = 0; i < order.size(); ++i) {
        c.addNode(Node{order[i], i == 0 ? -1 : index(parent.at(order[i])),
                       i == 0 ? kNoTable : node_table.at(order[i])});
    }

    // Wells of the right phase whose group is a node of this network.
    for (const auto& record : deck["WCONINJE"].front()) {
        if (record.getItem("TYPE").getTrimmedString(0) != well_phase) {
            continue;
        }
        const auto name = record.getItem("WELL").getTrimmedString(0);
        const auto group = well_group.at(name);
        if (!node_table.count(group)) {
            continue;
        }
        Well w;
        w.name = name;
        w.node = index(group);
        w.vfp_table = record.getItem("VFP_TABLE").get<int>(0);
        // BHP and RATE are UDA items even when the deck gives them as numbers.
        // RATE carries no usable dimension -- which unit it is in depends on the
        // injected phase, which the parser cannot know -- so getSI() hands back
        // the raw deck number and the conversion has to be done here. Getting
        // this wrong is silent: the limit comes out 86400x too large and the
        // rate control simply never activates.
        w.bhp_limit = record.getItem("BHP").get<UDAValue>(0).getSI();
        w.rate_limit = deck.getActiveUnitSystem().to_si(
            fluid == Fluid::Gas ? UnitSystem::measure::gas_surface_rate
                                : UnitSystem::measure::liquid_surface_rate,
            record.getItem("RATE").get<UDAValue>(0).get<double>());
        c.addWell(w);
    }

    for (const auto& keyword : deck["VFPINJ"]) {
        c.addInjTable(keyword);
    }
    return c;
}

Reference Reference::fromSummary(const std::string& prefix,
                                 const NetworkCase& c,
                                 const double time_days)
{
    const EclIO::ESmry summary(prefix + ".SMSPEC");
    // Report steps, not the raw vectors: a reference run's summary carries
    // ministeps too, and an operating point taken at one of those is a slightly
    // different state from the report the deck asked for.
    const auto time = summary.get_at_rstep("TIME");
    std::size_t at = 0;
    for (std::size_t i = 1; i < time.size(); ++i) {
        if (std::abs(time[i] - time_days) < std::abs(time[at] - time_days)) {
            at = i;
        }
    }

    const std::string rate_key = (c.fluid() == Fluid::Gas) ? "WGIR:" : "WWIR:";

    Reference ref;
    for (const auto& w : c.wells()) {
        ref.set(w.name, summary.get_at_rstep(rate_key + w.name)[at],
                summary.get_at_rstep("WBHP:" + w.name)[at]);
    }
    return ref;
}

void NetworkCase::calibrate(const Reference& reference)
{
    for (auto& w : wells_) {
        if (reference.has(w.name)) {
            w.q_ref = reference.rate(w.name);
            w.bhp_ref = reference.bhp(w.name);
        }
    }
    finish();
}

/// The reference at day 31 (rate and bhp per well), read off
/// opm-tests/eclref/e100reference/GNETINJE_GAS-01_ECL.
Reference referenceGnetinjeGasDay31()
{
    Reference r;
    r.set("G-3H", 486500.2, 295.3923);
    r.set("G-4H", 486530.8, 295.3244);
    r.set("F-1H", 276481.3, 295.0190);
    r.set("F-2H", 277082.5, 294.9374);
    return r;
}

/// GNETINJE_GAS-01 with its reference point set by hand, so the bench needs no
/// files. This is the same thing fromDeck() + Reference::fromSummary() produce
/// for that case, and setting the point by hand is also how you would pose a
/// case that no reference run covers.
NetworkCase gnetinjeGas()
{
    const auto sm3_day = cubic(meter) / day;
    const double bhp_limit = convert::from(425.0, bars);     // WCONINJE
    const double rate_limit = convert::from(1.0e6, sm3_day); // WCONINJE

    NetworkCase c;
    c.setFluid(Fluid::Gas);
    c.addTable(vfp_well);
    c.addTable(vfp_m5n);
    c.addTable(vfp_m5s);

    c.setTerminalPressure(convert::from(340.0, bars));       // GNETINJE PLAT-A
    c.addNode(Node{"PLAT-A", -1, kNoTable});
    c.addNode(Node{"M5S", 0, 3});
    c.addNode(Node{"M5N", 1, 2});
    c.addNode(Node{"G1", 1, kNoTable});
    c.addNode(Node{"F1", 2, kNoTable});

    for (const auto& [name, node] : std::initializer_list<std::pair<const char*, int>>{
             {"G-3H", 3}, {"G-4H", 3}, {"F-1H", 4}, {"F-2H", 4}}) {
        Well w;
        w.name = name;
        w.node = node;
        w.vfp_table = 1;
        w.bhp_limit = bhp_limit;
        w.rate_limit = rate_limit;
        c.addWell(w);
    }

    c.setStiffness(6.0e4);
    c.calibrate(referenceGnetinjeGasDay31());
    return c;
}

// ---------------------------------------------------------------------------
// Two formulations of the same case
//
// Eliminated: the unknowns are the pressures of the nodes that carry a table.
// Every rate is recovered from them by an inner solve, so the residual is cheap
// to state and awkward to differentiate -- the control clamps put kinks in it.
//
// Full: node pressures, branch rates, each well's (rate, bhp) and, when a group
// target is active, its multiplier. The eliminations become equations. Bigger,
// smooth within an active set, and the shape a reservoir/well system could
// absorb.
//
// Both expose size(), residual(x), start(p) and limitStep(), so the Newton
// below does not know which one it is solving.
// ---------------------------------------------------------------------------

using State = std::vector<double>;

/// Residual scaling: pressure equations in bar, rate equations in units of
/// kRateScale. Without this the two kinds of row differ by ~7 decades and no
/// single convergence tolerance means anything.
const double kPressureScale = convert::from(1.0, bars);
const double kRateScale = convert::from(1.0e4, cubic(meter) / day);

class EliminatedProblem
{
public:
    explicit EliminatedProblem(const NetworkCase& c) : case_(c) {}

    static constexpr const char* name = "eliminated";
    int size() const { return static_cast<int>(case_.solvedNodes().size()); }

    /// The fixed-point map: applied node pressures in, computed ones out.
    State G(const State& applied) const
    {
        const auto p = case_.nodePressures(applied);
        const auto q = case_.branchFlows(case_.rates(p));

        State out(size());
        const auto& solved = case_.solvedNodes();
        for (std::size_t i = 0; i < solved.size(); ++i) {
            const auto& node = case_.nodes()[solved[i]];
            out[i] = case_.tableBhp(node.vfp_table, p[node.parent], q[solved[i]]);
        }
        return out;
    }

    State residual(const State& x) const
    {
        const auto g = G(x);
        State r(size());
        for (int i = 0; i < size(); ++i) {
            r[i] = (g[i] - x[i]) / kPressureScale;
        }
        return r;
    }

    State start(const State& p) const { return p; }
    State pressures(const State& x) const { return x; }
    /// The eliminated form has no rate unknowns; recover them from the case.
    State wellRates(const State& x) const { return case_.rates(case_.nodePressures(x)); }
    double columnScale(const int) const { return kPressureScale; }
    State limitStep(const State&, const State& dx) const { return dx; }

private:
    const NetworkCase& case_;
};

/// The same case without the eliminations, solved by the library's
/// NetworkSolve::System -- the very code the simulator runs. The bench only
/// supplies the wells' inflow performance from its reference operating point,
/// where the simulator supplies it from the well Jacobian.
class FullProblem
{
public:
    explicit FullProblem(const NetworkCase& c)
        : system_(c.system()), solved_(c.solvedNodes()), terminal_(c.terminalPressure())
    {}

    static constexpr const char* name = "full";

    int size() const { return system_.size(); }
    State residual(const State& x) const { return system_.residual(x); }
    bool updateControls(const State& x) { return system_.updateControls(x); }
    double columnScale(const int i) const { return system_.columnScale(i); }
    State limitStep(const State& x, const State& dx) const
    {
        return enforce_bounds_ ? system_.limitStep(x, dx) : dx;
    }

    void setEnforceBounds(const bool on) { enforce_bounds_ = on; }
    void setAnalyticJacobian(const bool on) { system_.setAnalyticJacobian(on); }
    void setGuidesFromPotential(const bool on) { system_.setGuidesFromPotential(on); }
    void dropLastFromGroup() { system_.dropLastFromGroup(); }
    State wellRates(const State& x) const { return system_.wellRates(x); }
    const NetworkSolve::InjectionSystem<double>& system() const { return system_; }
    NetworkSolve::InjectionSystem<double>& system() { return system_; }

    /// The bench starts both formulations from the same applied node pressures.
    State start(const State& applied) const
    {
        return system_.start(applied_to_all_nodes_(applied));
    }

    /// Only the nodes the eliminated form solves for, so the two are comparable.
    State pressures(const State& x) const
    {
        const auto all = system_.pressures(x);
        State p;
        for (const int n : solved_) {
            p.push_back(all[n]);
        }
        return p;
    }

private:
    State applied_to_all_nodes_(const State& applied) const
    {
        State p(system_.nodes().size(), terminal_);
        for (std::size_t n = 1; n < system_.nodes().size(); ++n) {
            const auto it = std::find(solved_.begin(), solved_.end(), static_cast<int>(n));
            p[n] = (it != solved_.end()) ? applied[it - solved_.begin()]
                                         : p[system_.nodes()[n].parent];
        }
        return p;
    }

    NetworkSolve::InjectionSystem<double> system_;
    std::vector<int> solved_;
    double terminal_ = 0.0;
    bool enforce_bounds_ = false;
};

// ---------------------------------------------------------------------------
// Solvers
//
// Everything works on the residual F(x). Convergence is in the max norm of the
// scaled residual; the trust region uses the 2-norm because that is what its
// reduction ratio is defined against.
// ---------------------------------------------------------------------------

// Start where the simulator does: the wells' WCONINJE THP.
const State kStart{convert::from(400.0, bars), convert::from(400.0, bars)};
const double kTol = 0.01;                     // scaled, so 0.01 bar
const double kMaxStep = 100.0;                // column scales, so 100 bar
constexpr int kMaxIter = 200;

struct Result
{
    bool converged = false;
    int iterations = 0;
    State p{};
    State well_rate{};
};

double normMax(const State& v)
{
    double m = 0.0;
    for (const double e : v) {
        m = std::max(m, std::abs(e));
    }
    return m;
}
double norm2(const State& v)
{
    double sum = 0.0;
    for (const double e : v) {
        sum += e * e;
    }
    return std::sqrt(sum);
}
State operator+(const State& a, const State& b)
{
    State c(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        c[i] = a[i] + b[i];
    }
    return c;
}
State operator-(const State& a, const State& b)
{
    State c(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        c[i] = a[i] - b[i];
    }
    return c;
}
State operator*(const double a, const State& v)
{
    State c(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        c[i] = a * v[i];
    }
    return c;
}
State operator-(const State& v) { return -1.0 * v; }

// --- fixed-point methods, eliminated form only -------------------------------

Result damped(const EliminatedProblem& problem, State p, const double omega)
{
    for (int it = 1; it <= kMaxIter; ++it) {
        const auto r = problem.residual(p);
        if (normMax(r) < kTol) {
            return {true, it, p};
        }
        for (int i = 0; i < problem.size(); ++i) {
            p[i] = NodePressureUpdater<double>::damped(p[i], r[i] * kPressureScale, omega,
                                                       kMaxStep * kPressureScale);
        }
    }
    return {false, kMaxIter + 1, p};
}

Result bracketing(const EliminatedProblem& problem, State p, const double omega)
{
    std::vector<NodePressureUpdater<double>> updater(problem.size());
    for (int it = 1; it <= kMaxIter; ++it) {
        const auto g = problem.G(p);
        if (normMax(g - p) < kTol * kPressureScale) {
            return {true, it, p};
        }
        for (int i = 0; i < problem.size(); ++i) {
            p[i] = updater[i].next(p[i], g[i], /*valid=*/true, omega, kMaxStep * kPressureScale);
        }
    }
    return {false, kMaxIter + 1, p};
}

// --- Newton ------------------------------------------------------------------

/// Dense square system, small enough that Gaussian elimination with partial
/// pivoting is the whole story.
class Matrix
{
public:
    explicit Matrix(const int n) : n_(n), a_(n * n, 0.0) {}

    double& operator()(const int i, const int j) { return a_[i * n_ + j]; }

    /// Solves A y = b. Returns false if A is singular to working precision.
    bool solve(State b, State& y) const
    {
        auto a = a_;
        y.assign(n_, 0.0);
        for (int k = 0; k < n_; ++k) {
            int pivot = k;
            for (int i = k + 1; i < n_; ++i) {
                if (std::abs(a[i * n_ + k]) > std::abs(a[pivot * n_ + k])) {
                    pivot = i;
                }
            }
            if (std::abs(a[pivot * n_ + k]) < 1e-300) {
                return false;
            }
            if (pivot != k) {
                for (int j = 0; j < n_; ++j) {
                    std::swap(a[k * n_ + j], a[pivot * n_ + j]);
                }
                std::swap(b[k], b[pivot]);
            }
            for (int i = k + 1; i < n_; ++i) {
                const double f = a[i * n_ + k] / a[k * n_ + k];
                for (int j = k; j < n_; ++j) {
                    a[i * n_ + j] -= f * a[k * n_ + j];
                }
                b[i] -= f * b[k];
            }
        }
        for (int i = n_ - 1; i >= 0; --i) {
            double sum = b[i];
            for (int j = i + 1; j < n_; ++j) {
                sum -= a[i * n_ + j] * y[j];
            }
            y[i] = sum / a[i * n_ + i];
        }
        return true;
    }

private:
    int n_;
    std::vector<double> a_;
};

/// Finite-difference Jacobian. On the eliminated problem the well response has
/// kinks where a control switches, so a difference taken across one is not the
/// local slope -- which is exactly why the step needs globalising. The full
/// problem holds its controls fixed while this is taken, so it has no kinks.
template <class Problem>
Matrix jacobian(const Problem& problem, const State& x, const State& r)
{
    const int n = problem.size();
    Matrix J(n);
    for (int j = 0; j < n; ++j) {
        State shifted = x;
        const double h = 1e-2 * problem.columnScale(j);
        shifted[j] += h;
        const auto rj = problem.residual(shifted);
        for (int i = 0; i < n; ++i) {
            J(i, j) = (rj[i] - r[i]) / h;
        }
    }
    return J;
}

// --- globalisation strategies ------------------------------------------------
//
// Each takes the current point, the residual there and the full Newton step, and
// returns the point to move to. The Newton direction is the same in all of them.

/// Take the step as it comes.
struct FullStep
{
    static constexpr const char* name = "newton, full step";

    template <class Problem>
    State accept(const Problem&, const State& x, const State&, const State& dx)
    {
        return x + dx;
    }
};

/// Clamp each component, the way --network-max-pressure-update-in-bars does.
struct CappedStep
{
    static constexpr const char* name = "newton, capped step";
    double cap = kMaxStep;   // in column scales, so 100 means 100 bar

    template <class Problem>
    State accept(const Problem& problem, const State& x, const State&, const State& dx)
    {
        State capped = dx;
        for (int i = 0; i < problem.size(); ++i) {
            const double limit = cap * problem.columnScale(i);
            capped[i] = std::clamp(capped[i], -limit, limit);
        }
        return x + capped;
    }
};

/// Backtrack until the residual norm drops. With sufficient_decrease it is the
/// Armijo condition rather than plain decrease.
struct LineSearch
{
    static constexpr const char* name = "newton, line search";
    bool sufficient_decrease = false;
    int max_halvings = 12;

    template <class Problem>
    State accept(const Problem& problem, const State& x, const State& r, const State& dx)
    {
        const double f0 = norm2(r);
        double lambda = 1.0;
        for (int k = 0; k < max_halvings; ++k) {
            const State trial = x + lambda * dx;
            const double f = norm2(problem.residual(trial));
            const double target = sufficient_decrease ? (1.0 - 1e-4 * lambda) * f0 : f0;
            if (f < target) {
                return trial;
            }
            lambda *= 0.5;
        }
        return x + lambda * dx;
    }
};

/// Classic trust region on ||F||_2: shrink the step to the radius, accept on the
/// ratio of actual to predicted reduction, and shrink and retry when it is poor.
struct TrustRegion
{
    static constexpr const char* name = "newton, trust region";
    double radius = 50.0;          // in column scales, so 50 means 50 bar
    double radius_max = 400.0;
    double radius_min = 1e-4;
    double radius_start = 50.0;

    template <class Problem>
    State accept(const Problem& problem, const State& x, const State& r, const State& dx)
    {
        const double f0 = norm2(r);
        // Measure the step in column scales, so one radius covers pressures and rates.
        State scaled(dx.size());
        for (int i = 0; i < problem.size(); ++i) {
            scaled[i] = dx[i] / problem.columnScale(i);
        }
        const double len = norm2(scaled);

        while (radius > radius_min) {
            const double lambda = (len > radius && len > 0.0) ? radius / len : 1.0;
            const State trial = x + lambda * dx;
            const double f = norm2(problem.residual(trial));
            // The step solves J dx = -r, so the linear model predicts (1 - lambda)*r.
            const double predicted = lambda * f0;
            const double rho = predicted > 0.0 ? (f0 - f) / predicted : -1.0;

            if (rho > 0.1) {
                if (rho > 0.75 && lambda < 1.0) {
                    radius = std::min(2.0 * radius, radius_max);
                }
                return trial;
            }
            radius *= 0.5;
        }
        // The region collapsed, which here means the Jacobian was taken across a
        // control switch. Take the smallest step and reopen rather than stalling.
        const double lambda = std::min(1.0, radius_min / std::max(len, radius_min));
        radius = radius_start;
        return x + lambda * dx;
    }
};

/// Newton on either formulation. The full problem reselects its well controls
/// once per iteration; converging with a control still moving is not converged.
template <class Problem, class Globalisation>
Result newton(Problem problem, const State& start, Globalisation g = {})
{
    State x = problem.start(start);
    for (int it = 1; it <= kMaxIter; ++it) {
        bool controls_moved = false;
        if constexpr (requires { problem.updateControls(x); }) {
            controls_moved = problem.updateControls(x);
        }
        const auto r = problem.residual(x);
        if (normMax(r) < kTol && !controls_moved) {
            return {true, it, problem.pressures(x), problem.wellRates(x)};
        }
        State dx;
        if (!jacobian(problem, x, r).solve(-r, dx)) {
            return {false, kMaxIter + 1, problem.pressures(x)};
        }
        // Keep the iterate inside the box the tables describe before anything
        // else looks at the step.
        dx = problem.limitStep(x, dx);
        // The residual jumps when a control switches, and that jump is not a
        // failure to make progress. Letting a globalisation veto it stalls the
        // active set instead of resolving it.
        x = controls_moved ? x + dx : g.accept(problem, x, r, dx);
    }
    return {false, kMaxIter + 1, problem.pressures(x)};
}

} // anonymous namespace



namespace {
    const State kExpected{convert::from(209.30, bars), convert::from(204.19, bars)};

    void report(const char* name, const Result& r)
    {
        BOOST_TEST_MESSAGE(std::left << std::setw(26) << name
                           << (r.converged ? "converged in " : "FAILED after ")
                           << std::setw(4) << r.iterations << " iterations, p = ("
                           << convert::to(r.p[0], bars) << ", "
                           << convert::to(r.p[1], bars) << ") bar");
    }

    /// The grid the basin tests sweep: starting node pressures across the tables' THP axis.
    std::vector<State> startingPoints()
    {
        std::vector<State> starts;
        for (int a = 60; a <= 500; a += 20) {
            for (int b = 60; b <= 500; b += 20) {
                starts.push_back({convert::from(a, bars), convert::from(b, bars)});
            }
        }
        return starts;
    }

    /// How many of the starting points a method reaches the right answer from.
    template <class Solve>
    int basin(const char* name, Solve&& solve, const State& expected = kExpected)
    {
        const auto starts = startingPoints();
        const double tol = convert::from(1.0, bars);

        int solved = 0, iterations = 0;
        for (const auto& start : starts) {
            const auto r = solve(start);
            if (r.converged && normMax(r.p - expected) < tol) {
                ++solved;
                iterations += r.iterations;
            }
        }
        BOOST_TEST_MESSAGE(std::left << std::setw(26) << name << solved << "/" << starts.size()
                           << " starts, mean " << (solved ? iterations / solved : 0)
                           << " iterations");
        return solved;
    }
}



namespace {
    /// Where opm-tests lives, if it does. The reference-driven cases skip
    /// without it, since it is not a build dependency.
    const std::string kTests = "/Users/hnil/Documents/OPM/opm_feature/opm-tests/";

    struct DeckCase
    {
        std::string deck;
        std::string reference;
        Fluid fluid;
        std::string rate_key;     // field injection rate
    };

    const DeckCase kGasCase{kTests + "network/GNETINJE_GAS-01.DATA",
                            kTests + "eclref/e100reference/GNETINJE_GAS-01_ECL",
                            Fluid::Gas, "FGIR"};
    const DeckCase kWaterCase{kTests + "network/GNETINJE_WAT-01.DATA",
                              kTests + "eclref/e100reference/GNETINJE_WAT-01_ECL",
                              Fluid::Water, "FWIR"};

    bool available(const DeckCase& c)
    {
        return std::filesystem::exists(c.deck) && std::filesystem::exists(c.reference + ".SMSPEC");
    }

    /// Solve the case at every report step of its reference run and report how
    /// far the node pressures land from it. Returns (matched, considered).
    ///
    /// At each step the wells are calibrated to what the reference says they
    /// were doing, and the field target is applied when the reference has the
    /// wells on group control -- WMCTL 6 is THP, anything negative is a group.
    /// So this asks whether the network side reproduces the reference at operating
    /// points across the whole run, not just the one the bench was built on.
    std::pair<int, int> sweepReference(const DeckCase& deck_case, const double tol_bar)
    {
        const EclIO::ESmry summary(deck_case.reference + ".SMSPEC");
        const auto time = summary.get_at_rstep("TIME");
        const auto field_rate = summary.get_at_rstep(deck_case.rate_key);
        const bool gas = deck_case.fluid == Fluid::Gas;
        const auto sm3d = cubic(meter) / day;

        const auto base = fromDeck(deck_case.deck, deck_case.fluid);
        const auto control = summary.get_at_rstep("WMCTL:" + base.wells().front().name);

        int matched = 0, considered = 0;
        for (std::size_t step = 0; step < time.size(); ++step) {
            if (time[step] <= 0.0 || field_rate[step] <= 0.0) {
                continue;      // nothing injected yet
            }
            ++considered;

            auto c = fromDeck(deck_case.deck, deck_case.fluid);
            c.calibrate(Reference::fromSummary(deck_case.reference, c, time[step]));
            c.setRelativeStiffness(0.12);
            const bool on_group = control[step] < 0.0;
            if (on_group) {
                c.setGroupTarget(convert::from(field_rate[step], sm3d));
            }
            c.finish();

            const auto solved = newton(FullProblem{c}, kStart, FullStep{});
            std::string detail;
            double worst = 0.0;
            for (std::size_t i = 0; i < c.solvedNodes().size(); ++i) {
                const auto& node = c.nodes()[c.solvedNodes()[i]].name;
                const double reference =
                    summary.get_at_rstep((gas ? "GPRG:" : "GPRW:") + node)[step];
                const double got = solved.converged ? convert::to(solved.p[i], bars) : 0.0;
                worst = std::max(worst, std::abs(got - reference));
                detail += fmt::format(" {} {:.2f}/{:.2f}", node, got, reference);
            }
            const bool ok = solved.converged && worst < tol_bar;
            matched += ok ? 1 : 0;
            BOOST_TEST_MESSAGE(fmt::format("  t={:6.1f} {:5} {:11} {:<26} worst {:.2f} bar {}",
                                           time[step], on_group ? "GRUP" : "THP",
                                           solved.converged
                                               ? fmt::format("{} it", solved.iterations) : "FAILED",
                                           detail, worst, ok ? "" : "  <--"));
        }
        BOOST_TEST_MESSAGE("  " << matched << "/" << considered << " report steps within "
                           << tol_bar << " bar");
        return {matched, considered};
    }
}




// The prototype's network, built to order so several cases can share it.
// FIELD -> PROD at 80 bar, two producers of different productivity on the node.
// The table has to outlive the properties object, which keeps a reference.
class ProductionCase
{
public:
    using Sys = NetworkSolve::ProductionSystem<double>;

    explicit ProductionCase(const double productivity_2 = 0.7)
    {
        const auto deck = Parser{}.parseString(vfp_prod);
        tables_.emplace_back(deck["VFPPROD"].front(), /*gaslift_opt_active=*/false, UnitSystem{});
        props_.addTable(tables_.back());

        for (const auto& [name, productivity] :
             std::initializer_list<std::pair<const char*, double>>{{"P-1", 1.0},
                                                                   {"P-2", productivity_2}}) {
            Sys::Well w;
            w.name = name;
            w.node = 1;
            w.vfp_table = 3;
            w.bhp_limit = convert::from(40.0, bars);
            w.oil_rate_limit = convert::from(600.0, cubic(meter) / day);
            w.in_group = true;
            const double q0 = convert::from(400.0 * productivity, cubic(meter) / day);
            const double slope = q0 / convert::from(120.0, bars);
            for (int ph = 0; ph < Sys::NP; ++ph) {
                const double share = (ph == 0) ? 0.3 : (ph == 1) ? 0.7 : 70.0;
                w.ipr_a[ph] = share * q0 * 2.0;
                w.ipr_b[ph] = -share * slope * 2.0;
            }
            wells_.push_back(w);
        }
    }

    std::vector<Sys::Well>& wells() { return wells_; }
    void setGroupTarget(const double t) { target_ = t; }
    double groupTarget() const { return target_; }

    Sys system() const
    {
        Sys s(props_, units_);
        s.setTerminalPressure(convert::from(80.0, bars));
        s.addNode(NetworkSolve::Node{"FIELD", -1, NetworkSolve::NoTable}, 0.0);
        s.addNode(NetworkSolve::Node{"PROD", 0, 3}, 0.0);
        for (const auto& w : wells_) {
            s.addWell(w);
        }
        s.setGroupTarget(target_);
        s.finish();
        return s;
    }

    /// What the two wells produce with no group target, at the pressures they
    /// settle on themselves. The scale every target here is quoted against.
    double freeTotal() const
    {
        ProductionCase open(*this);
        open.target_ = 0.0;
        auto s = open.system();
        const auto r = NetworkSolve::solve(s, guess(), kParams, NetworkSolve::FullStep{});
        BOOST_REQUIRE(r.converged);
        return r.well_rate[0] + r.well_rate[1];
    }

    static std::vector<double> guess()
    {
        return {convert::from(80.0, bars), convert::from(90.0, bars)};
    }

private:
    std::deque<VFPProdTable> tables_;
    VFPProdProperties<double> props_;
    UnitSystem units_{};
    std::vector<Sys::Well> wells_;
    double target_ = 0.0;
};



// The group tree as sparse rows, two levels deep.
//
//        FIELD (target)
//        /          \
//      P1            P2
//    W1  W2          W3
//
// Everything measured so far had exactly one group constraint, flattened
// outside the Newton into one target and per-well guide rates. This writes the
// tree out: phase rates and a multiplier per group, a definition row that
// touches only immediate children, and one complementarity per group saying it
// sits at whichever of its parent's share and its own target is smaller.
//
// The well rows are four FB deep here (rate, tubing, bhp, share), which is the
// nesting depth nothing has tested.
namespace {
    struct TreeCase
    {
        using Sys = NetworkSolve::ProductionSystem<double>;
        VFPProdProperties<double> props;
        UnitSystem units{};

        // bhp limit 100 bar, shut-in 250 bar, cap = rate at the bhp limit.
        static std::array<double, 2> ipr(const double cap_per_day)
        {
            const double q = cap_per_day / 86400.0;
            const double b = -q / (250e5 - 100e5);
            return {-b * 250e5, b};
        }

        /// water_cut: every well makes this fraction of its oil as water, so
        /// a liquid target of (1 + wct) * T must give the same answer as an
        /// oil target of T. That equivalence is the multi-phase check.
        Sys build(const double field_target, const double p1_target,
                  const std::array<double, 3>& caps,
                  const typename Sys::Mode mode = Sys::Mode::Oil,
                  const double water_cut = 0.0,
                  const bool active_set = false) const
        {
            Sys sys(props, units);
            NetworkSolve::Node root; root.name = "TERM"; root.parent = -1;
            sys.addNode(root, 0.0);
            NetworkSolve::Node n; n.name = "N"; n.parent = 0;
            sys.addNode(n, 0.0);
            sys.setTerminalPressure(convert::from(120.0, bars));

            typename Sys::Group f; f.name = "FIELD"; f.parent = -1;
            f.target = field_target / 86400.0;
            f.mode = mode;
            const int gf = sys.addGroup(f);
            typename Sys::Group p1; p1.name = "P1"; p1.parent = gf; p1.guide = 1.0;
            if (p1_target > 0.0) { p1.target = p1_target / 86400.0; }
            const int g1 = sys.addGroup(p1);
            typename Sys::Group p2; p2.name = "P2"; p2.parent = gf; p2.guide = 1.0;
            const int g2 = sys.addGroup(p2);

            const int grp[3] = {g1, g1, g2};
            for (int i = 0; i < 3; ++i) {
                typename Sys::Well w;
                w.name = "W" + std::to_string(i + 1);
                w.node = 1;
                w.vfp_table = 0;                 // no tubing: rate/bhp/share only
                const auto ab = ipr(caps[i]);
                w.ipr_a[1] = ab[0]; w.ipr_b[1] = ab[1];
                // Water on the same line, so the cut is constant in bhp.
                w.ipr_a[0] = water_cut * ab[0]; w.ipr_b[0] = water_cut * ab[1];
                w.bhp_limit = convert::from(100.0, bars);
                w.oil_rate_limit = 5000.0 / 86400.0;
                w.guide = 1.0;
                w.group = grp[i];
                w.q_start = 0.5 * caps[i] / 86400.0;
                sys.addWell(std::move(w));
            }
            sys.setGroupTree(true);
            sys.setAnalyticJacobian(true);
            sys.setComplementarity(true);
            sys.setGroupActiveSet(active_set);
            sys.finish();
            sys.finishGroups();
            return sys;
        }
    };
}



// ---------------------------------------------------------------------------
// The group tree on a real network.
//
// Everything above solved the tree with no tubing, so the node pressure moved
// nothing. Here the tree of GROUPTREE.DATA -- FIELD -> PLAT -> {GP1 -> W1 W2,
// GP2 -> W3} -- hangs on branches with the VFPPROD table, so a well's capacity
// depends on the node pressure, the node pressure on the allocation, and the
// two have to be found together. Three routes to the same answer:
// complementarity inside the Newton, the tree pass every iteration, and the
// tree pass between converged solves (NetworkTreeSolve.hpp).
// ---------------------------------------------------------------------------
namespace {
    struct NetTreeCase
    {
        using Sys = NetworkSolve::ProductionSystem<double>;
        // props keeps a reference to the table, so the table lives here too.
        VFPProdTable table;
        VFPProdProperties<double> props;
        UnitSystem units{};
        NetTreeCase()
            : table(Parser{}.parseString(vfp_prod)["VFPPROD"].front(), false, UnitSystem{})
        {
            props.addTable(table);
        }

        /// Liquid q(bhp) = 2 q0 (1 - bhp / shut_in), split 0.3 water, 0.7 oil,
        /// GOR 100 -- the prototype's wells. water_shut_in moves the water
        /// line's zero away from the oil line's, so the cut varies with bhp.
        static void ipr(Sys::Well& w, const double q0_sm3d, const double water_shut_in_bar)
        {
            const double q0 = q0_sm3d / 86400.0;
            const double oil_shut = convert::from(120.0, bars);
            const double water_shut = convert::from(water_shut_in_bar, bars);
            const double share[3] = {0.3, 0.7, 70.0};   // water, oil, gas
            for (int ph = 0; ph < Sys::NP; ++ph) {
                const double shut = (ph == 0) ? water_shut : oil_shut;
                w.ipr_a[ph] = 2.0 * share[ph] * q0;
                w.ipr_b[ph] = -2.0 * share[ph] * q0 / shut;
            }
        }

        /// fork: N1 and N2 both hang off the terminal; chain: N2 under N1.
        Sys build(const double plat_target, const std::array<double, 3>& q0,
                  const bool fork, const Sys::Mode mode = Sys::Mode::Oil,
                  const double water_shut_in_bar = 120.0) const
        {
            Sys sys(props, units);
            sys.addNode(NetworkSolve::Node{"TERM", -1, NetworkSolve::NoTable}, 0.0);
            sys.addNode(NetworkSolve::Node{"N1", 0, 3}, 0.0);
            sys.addNode(NetworkSolve::Node{"N2", fork ? 0 : 1, 3}, 0.0);
            sys.setTerminalPressure(convert::from(20.0, bars));

            typename Sys::Group f; f.name = "FIELD"; f.parent = -1;
            const int gf = sys.addGroup(f);
            typename Sys::Group plat; plat.name = "PLAT"; plat.parent = gf;
            plat.target = plat_target / 86400.0; plat.mode = mode; plat.guide = 3.0;
            const int gp = sys.addGroup(plat);
            typename Sys::Group g1; g1.name = "GP1"; g1.parent = gp; g1.guide = 2.0;
            const int gp1 = sys.addGroup(g1);
            typename Sys::Group g2; g2.name = "GP2"; g2.parent = gp; g2.guide = 1.0;
            const int gp2 = sys.addGroup(g2);

            const int grp[3] = {gp1, gp1, gp2};
            const int node[3] = {1, 1, 2};
            for (int i = 0; i < 3; ++i) {
                typename Sys::Well w;
                w.name = "W" + std::to_string(i + 1);
                w.node = node[i];
                w.vfp_table = 3;
                ipr(w, q0[i], water_shut_in_bar);
                w.bhp_limit = convert::from(40.0, bars);
                w.oil_rate_limit = 5000.0 / 86400.0;
                w.guide = 1.0;
                w.group = grp[i];
                w.q_start = 0.5 * std::max(w.ipr_a[1] + w.ipr_b[1] * w.bhp_limit, 0.0);
                sys.addWell(std::move(w));
            }
            sys.setGroupTree(true);
            sys.setAnalyticJacobian(true);
            sys.setComplementarity(true);
            sys.finish();
            sys.finishGroups();
            return sys;
        }

        static std::vector<double> guess(const double bar)
        {
            return {convert::from(20.0, bars), convert::from(bar, bars), convert::from(bar, bars)};
        }
    };

    /// Stein's balancer fed with the system's own capacities at a state: each
    /// well's limit is what its own controls allow at the node pressure there,
    /// its phase rates the IPR at the bhp that allowance implies. Same IPR,
    /// same tubing curve, same pressures as the Newton. The tree, the targets
    /// and the guides are the system's; the Schedule only supplies GuideRate.
    /// Returns oil rates per well in sm3/d, or empty if the balancer rejected
    /// the tree.
    std::vector<double> steinAllocation(NetworkSolve::ProductionSystem<double>& sys,
                                        const std::vector<double>& x,
                                        const Schedule& schedule, const int step = 0)
    {
        using Sys = NetworkSolve::ProductionSystem<double>;
        sys.updateControls(x);      // the capacities at this state's pressures
        GuideRate guide_rate{schedule};
        DeferredLogger logger;
        // Potentials for every well and, summed beneath, every group: a
        // guide rate on a mode other than oil -- a liquid target over oil
        // children -- is read off the potentials, and a group with none is
        // a lookup failure inside the balancer.
        const auto& groups = sys.groups();
        std::vector<std::array<double, 3>> gpot(sys.numGroups(), {0.0, 0.0, 0.0});   // oil, gas, water
        for (int w = 0; w < sys.numWells(); ++w) {
            // The well's guide as its oil potential, the other phases in the
            // well's proportions at its bhp limit, so a group's guide on any
            // mode is the sum of its wells' -- the system's own rule.
            const auto& well = sys.wells()[w];
            const double bhp = well.bhp_limit;
            const double oil = std::max(well.ipr_a[1] + well.ipr_b[1] * bhp, 1e-30);
            // A well the dump gave no guide (individually controlled) still
            // needs one the balancer can divide by.
            const double guide = std::max(well.guide, 1e-9 / 86400.0);
            const std::array<double, 3> pot{guide,
                                            guide * std::max(well.ipr_a[2] + well.ipr_b[2] * bhp, 0.0) / oil,
                                            guide * std::max(well.ipr_a[0] + well.ipr_b[0] * bhp, 0.0) / oil};
            guide_rate.compute(well.name, step, 0.0, pot[0], pot[1], pot[2]);
            for (int g = well.group; g >= 0; g = groups[g].parent) {
                for (int i = 0; i < 3; ++i) { gpot[g][i] += well.efficiency * pot[i]; }
            }
        }
        for (int g = 0; g < sys.numGroups(); ++g) {
            guide_rate.compute(groups[g].name, step, 0.0, gpot[g][0], gpot[g][1], gpot[g][2]);
        }
        auto wmode = [](const Sys::Mode m) {
            switch (m) {
            case Sys::Mode::Water:  return Opm::Well::ProducerCMode::WRAT;
            case Sys::Mode::Gas:    return Opm::Well::ProducerCMode::GRAT;
            case Sys::Mode::Liquid: return Opm::Well::ProducerCMode::LRAT;
            default:                return Opm::Well::ProducerCMode::ORAT;
            }
        };
        auto gmode = [](const Sys::Mode m) {
            switch (m) {
            case Sys::Mode::Water:  return Opm::Group::ProductionCMode::WRAT;
            case Sys::Mode::Gas:    return Opm::Group::ProductionCMode::GRAT;
            case Sys::Mode::Liquid: return Opm::Group::ProductionCMode::LRAT;
            default:                return Opm::Group::ProductionCMode::ORAT;
            }
        };
        ProdGroupTreeBalancer::Tree<double> tree;
        for (int g = 0; g < sys.numGroups(); ++g) {
            ProdGroupTreeNode<double> n;
            n.name = groups[g].name;
            n.type = ProdNodeType::Group;
            n.parent = groups[g].parent >= 0 ? groups[groups[g].parent].name : std::string{};
            n.availableForGroupControl = true;
            n.efficiencyFactor = groups[g].efficiency;
            // As populateGroupNode(): a group has a guide rate only if the
            // deck defines one (GCONPROD item 9); none of these do.
            n.hasGuideRate = false;
            const bool has_target = groups[g].target > 0.0;
            n.modeCategory = has_target ? ProdNodeModeCategory::Individual
                                        : ProdNodeModeCategory::Group;
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
            ProdGroupTreeNode<double> n;
            n.name = well.name;
            n.type = ProdNodeType::Well;
            n.parent = groups[well.group].name;
            n.availableForGroupControl = true;
            n.efficiencyFactor = well.efficiency;
            n.mode = Opm::Well::ProducerCMode::GRUP;
            n.hasGuideRate = true;
            const double allow = sys.ownAllowance(w);
            const double bhp = (allow - well.ipr_a[1]) / well.ipr_b[1];
            n.Limits[Opm::Well::ProducerCMode::ORAT] = allow;
            // canonical [oil, water, gas], production negative
            n.rates = {-std::max(well.ipr_a[1] + well.ipr_b[1] * bhp, 0.0),
                       -std::max(well.ipr_a[0] + well.ipr_b[0] * bhp, 0.0),
                       -std::max(well.ipr_a[2] + well.ipr_b[2] * bhp, 0.0)};
            n.initialRates = n.rates;
            tree.at(n.parent).children.push_back(n.name);
            tree.emplace(well.name, std::move(n));
        }
        if (!ProdGroupTreeBalancer::balanceTreeForTesting(tree, guide_rate, 1e-8, logger)) {
            return {};
        }
        std::vector<double> q(sys.numWells());
        for (int w = 0; w < sys.numWells(); ++w) {
            q[w] = -tree.at(sys.wells()[w].name).rates[0] * 86400.0;
        }
        return q;
    }

    std::string rateList(const NetTreeCase::Sys& sys, const std::vector<double>& q_sm3s)
    {
        std::string s;
        for (int w = 0; w < sys.numWells(); ++w) {
            s += fmt::format(" {}={:.1f}", sys.wells()[w].name, q_sm3s[w] * 86400.0);
        }
        return s;
    }
}



// ---------------------------------------------------------------------------
// Both trees from a real deck.
//
// The group tree (GRUPTREE / GCONPROD / WCONPROD) and the flow network
// (NODEPROP / BRANPROP / VFPPROD) of the same deck, as one system. They share
// only the wells: on MODEL5 the group C1 reports to M5N, which is no node,
// while its wells flow to PLAT-A directly. The IPR is not in the deck --
// there is no reservoir here -- so every producer gets a straight line,
// liquid J (p_res - bhp) split by a water cut and a GOR, and a well's other
// rate limits are turned into the oil rate at which that line meets them,
// which is exact for a linear IPR.
// ---------------------------------------------------------------------------
namespace {
    struct DeckTrees
    {
        using Sys = NetworkSolve::ProductionSystem<double>;
        struct Ipr {
            double p_res_bar = 280.0;
            /// oil at the bhp limit, as a multiple of the well's own ORAT limit
            double j_scale = 2.0;
            double wct = 0.2;
            double gor = 100.0;
            /// per-well multiplier on j_scale, for wells meant to be weak
            std::map<std::string, double> j_of;
        };
        std::unique_ptr<Deck> deck;
        std::unique_ptr<EclipseState> es;
        std::unique_ptr<Schedule> schedule;
        SummaryState st{TimeService::now(), 0.0};
        std::deque<VFPProdTable> tables;       // props keeps references
        VFPProdProperties<double> props;
        UnitSystem units{};
        std::vector<std::string> node_order;
        std::map<std::string, int> gidx;
        std::string description;

        struct FromText {};
        explicit DeckTrees(const std::string& path)
            : DeckTrees(Parser{}.parseFile(path))
        {}
        DeckTrees(const std::string& text, FromText)
            : DeckTrees(Parser{}.parseString(text))
        {}
        explicit DeckTrees(Deck parsed)
        {
            deck = std::make_unique<Deck>(std::move(parsed));
            es = std::make_unique<EclipseState>(*deck);
            schedule = std::make_unique<Schedule>(*deck, *es);
            for (const auto& kw : deck->getKeywordList("VFPPROD")) {
                tables.emplace_back(*kw, /*gaslift_opt_active=*/false, units);
                props.addTable(tables.back());
            }
        }

        Sys build(const int step, const Ipr& ipr)
        {
            const auto& sched = *schedule;
            const auto& network = sched[step].network();
            Sys sys(props, units);
            description.clear();

            // The flow network: the fixed-pressure root and everything under it.
            const auto roots = network.roots();
            BOOST_REQUIRE(!roots.empty());
            const auto& root = roots.front().get();
            BOOST_REQUIRE(root.terminal_pressure().has_value());
            sys.setTerminalPressure(*root.terminal_pressure());
            node_order = {root.name()};
            std::map<std::string, int> nidx{{root.name(), 0}};
            sys.addNode(NetworkSolve::Node{root.name(), -1, NetworkSolve::NoTable}, 0.0);
            for (std::size_t at = 0; at < node_order.size(); ++at) {
                for (const auto& br : network.downtree_branches(node_order[at])) {
                    const auto& child = br.downtree_node();
                    const int table = br.vfp_table().has_value() && props.hasTable(*br.vfp_table())
                        ? *br.vfp_table() : NetworkSolve::NoTable;
                    nidx[child] = static_cast<int>(node_order.size());
                    node_order.push_back(child);
                    sys.addNode(NetworkSolve::Node{child, static_cast<int>(at), table}, 0.0);
                    description += fmt::format(" {}->{}{}", child, node_order[at],
                                               table == NetworkSolve::NoTable ? "" : fmt::format("[{}]", table));
                }
            }

            // Producers: which are open, what they may make, where they flow.
            struct W { std::string name, group; int node; Opm::Well::ProductionControls ctl;
                       double eff; int vfp; };
            std::vector<W> wells;
            for (const auto& wname : sched.wellNames(step)) {
                const auto& well = sched.getWell(wname, step);
                if (!well.isProducer() || well.getStatus() != Opm::Well::Status::OPEN) { continue; }
                // A well flows into the node named as its group, or the nearest
                // group above it that is a node.
                std::string g = well.groupName();
                while (!g.empty() && !nidx.count(g)) {
                    g = (g == "FIELD") ? std::string{} : sched.getGroup(g, step).parent();
                }
                if (g.empty()) { continue; }         // not in the network
                wells.push_back({wname, well.groupName(), nidx.at(g),
                                 well.productionControls(st), well.getEfficiencyFactor(),
                                 well.vfp_table_number()});
            }

            // Every well's line, and the oil it makes at its bhp limit -- its
            // guide, as a potential would be. Group guides are the sums beneath.
            std::map<std::string, double> wguide;
            auto line = [&](const W& w, Sys::Well& out) {
                const double p_res = convert::from(ipr.p_res_bar, bars);
                const double bhp_lim = w.ctl.bhp_limit > 0.0 ? w.ctl.bhp_limit : unit::barsa;
                const double orat = w.ctl.hasControl(WellProducerCMode::ORAT) && w.ctl.oil_rate > 0.0
                    ? w.ctl.oil_rate : 1000.0 / 86400.0;
                const double jm = ipr.j_of.count(w.name) ? ipr.j_of.at(w.name) : 1.0;
                const double j_oil = jm * ipr.j_scale * orat / (p_res - bhp_lim);
                const double j_liq = j_oil / (1.0 - ipr.wct);
                const double j[3] = {ipr.wct * j_liq, j_oil, ipr.gor * j_oil};   // water, oil, gas
                for (int ph = 0; ph < Sys::NP; ++ph) { out.ipr_a[ph] = j[ph] * p_res; out.ipr_b[ph] = -j[ph]; }
                out.bhp_limit = bhp_lim;
                // Each rate limit as the oil rate at which the line meets it.
                double allow = 0.0;
                auto limit = [&](const bool has, const double value, const Sys::Mode m) {
                    if (!has || !(value > 0.0)) { return; }
                    const auto c = Sys::modeWeights(m, {});
                    double ca = 0.0, cb = 0.0;
                    for (int ph = 0; ph < Sys::NP; ++ph) { ca += c[ph] * out.ipr_a[ph]; cb += c[ph] * out.ipr_b[ph]; }
                    const double bhp = (value - ca) / cb;
                    const double oil = out.ipr_a[1] + out.ipr_b[1] * bhp;
                    allow = (allow > 0.0) ? std::min(allow, oil) : oil;
                };
                limit(w.ctl.hasControl(WellProducerCMode::ORAT), w.ctl.oil_rate, Sys::Mode::Oil);
                limit(w.ctl.hasControl(WellProducerCMode::WRAT), w.ctl.water_rate, Sys::Mode::Water);
                limit(w.ctl.hasControl(WellProducerCMode::GRAT), w.ctl.gas_rate, Sys::Mode::Gas);
                limit(w.ctl.hasControl(WellProducerCMode::LRAT), w.ctl.liquid_rate, Sys::Mode::Liquid);
                out.oil_rate_limit = allow;
                return out.ipr_a[1] + out.ipr_b[1] * bhp_lim;
            };
            std::vector<Sys::Well> built(wells.size());
            for (std::size_t i = 0; i < wells.size(); ++i) {
                built[i].name = wells[i].name;
                built[i].node = wells[i].node;
                built[i].vfp_table = (wells[i].vfp > 0 && props.hasTable(wells[i].vfp)) ? wells[i].vfp : 0;
                built[i].efficiency = wells[i].eff;
                wguide[wells[i].name] = line(wells[i], built[i]);
                built[i].guide = wguide[wells[i].name];
                built[i].q_start = 0.5 * wguide[wells[i].name];
            }
            std::function<double(const std::string&)> subtreeGuide = [&](const std::string& name) {
                const auto& grp = sched.getGroup(name, step);
                double g = 0.0;
                for (const auto& child : grp.groups()) { g += subtreeGuide(child); }
                for (const auto& w : grp.wells()) { if (wguide.count(w)) { g += wguide.at(w); } }
                return g;
            };

            // The group tree, parents first, with each group's own target.
            gidx.clear();
            std::function<void(const std::string&, int)> addTree = [&](const std::string& name, const int parent) {
                const auto& grp = sched.getGroup(name, step);
                typename Sys::Group g;
                g.name = name;
                g.parent = parent;
                g.efficiency = grp.getGroupEfficiencyFactor();
                g.guide = subtreeGuide(name);
                if (grp.isProductionGroup()) {
                    const auto ctl = grp.productionControls(st);
                    using C = Opm::Group::ProductionCMode;
                    if (ctl.cmode == C::ORAT)      { g.mode = Sys::Mode::Oil;    g.target = ctl.oil_target; }
                    else if (ctl.cmode == C::LRAT) { g.mode = Sys::Mode::Liquid; g.target = ctl.liquid_target; }
                    else if (ctl.cmode == C::GRAT) { g.mode = Sys::Mode::Gas;    g.target = ctl.gas_target; }
                    else if (ctl.cmode == C::WRAT) { g.mode = Sys::Mode::Water;  g.target = ctl.water_target; }
                }
                const int me = sys.addGroup(std::move(g));
                gidx[name] = me;
                if (sys.groups()[me].target > 0.0) {
                    description += fmt::format(" {}:{:.0f}", name, sys.groups()[me].target * 86400.0);
                }
                for (const auto& child : grp.groups()) { addTree(child, me); }
            };
            addTree("FIELD", -1);
            for (std::size_t i = 0; i < wells.size(); ++i) {
                built[i].group = gidx.at(wells[i].group);
                description += fmt::format(" {}@{}/{}", wells[i].name, node_order[wells[i].node], wells[i].group);
                sys.addWell(std::move(built[i]));
            }
            sys.setGroupTree(true);
            sys.setAnalyticJacobian(true);
            sys.setComplementarity(true);
            sys.finish();
            sys.finishGroups();
            return sys;
        }

        std::vector<double> guess(const double bar) const
        {
            std::vector<double> p(node_order.size(), convert::from(bar, bars));
            return p;
        }

        /// The deck's group tree onto a system whose wells came from elsewhere
        /// -- a simulator dump, with its own IPRs and guides. Wells not in the
        /// schedule stay outside the tree.
        void attachTree(Sys& sys, const int step)
        {
            const auto& sched = *schedule;
            std::map<std::string, double> wguide;
            for (int w = 0; w < sys.numWells(); ++w) { wguide[sys.wells()[w].name] = sys.wells()[w].guide; }
            std::function<double(const std::string&)> subtreeGuide = [&](const std::string& name) {
                const auto& grp = sched.getGroup(name, step);
                double g = 0.0;
                for (const auto& child : grp.groups()) { g += subtreeGuide(child); }
                for (const auto& w : grp.wells()) { if (wguide.count(w)) { g += wguide.at(w); } }
                return g;
            };
            gidx.clear();
            description.clear();
            std::function<void(const std::string&, int)> addTree = [&](const std::string& name, const int parent) {
                const auto& grp = sched.getGroup(name, step);
                typename Sys::Group g;
                g.name = name;
                g.parent = parent;
                g.efficiency = grp.getGroupEfficiencyFactor();
                g.guide = subtreeGuide(name);
                if (grp.isProductionGroup()) {
                    const auto ctl = grp.productionControls(st);
                    using C = Opm::Group::ProductionCMode;
                    if (ctl.cmode == C::ORAT)      { g.mode = Sys::Mode::Oil;    g.target = ctl.oil_target; }
                    else if (ctl.cmode == C::LRAT) { g.mode = Sys::Mode::Liquid; g.target = ctl.liquid_target; }
                    else if (ctl.cmode == C::GRAT) { g.mode = Sys::Mode::Gas;    g.target = ctl.gas_target; }
                    else if (ctl.cmode == C::WRAT) { g.mode = Sys::Mode::Water;  g.target = ctl.water_target; }
                }
                const int me = sys.addGroup(std::move(g));
                gidx[name] = me;
                if (sys.groups()[me].target > 0.0) {
                    description += fmt::format(" {}:{:.0f}", name, sys.groups()[me].target * 86400.0);
                }
                for (const auto& child : grp.groups()) { addTree(child, me); }
            };
            addTree("FIELD", -1);
            for (int w = 0; w < sys.numWells(); ++w) {
                const auto& name = sys.wells()[w].name;
                if (!sched.hasWell(name, step)) { continue; }
                const auto it = gidx.find(sched.getWell(name, step).groupName());
                if (it != gidx.end()) { sys.setWellGroup(w, it->second); }
            }
            sys.setGroupTree(true);
            sys.finishGroups();
        }
    };

    const std::string kNetworkDecks = "/Users/hnil/Documents/OPM/opm_feature/opm-tests/network/";

    /// The three routes and the oracle on one built system; a line of report.
    void threeRoutesAndOracle(DeckTrees& dt, const DeckTrees::Ipr& ipr, const int step,
                              const double guess_bar, const char* tag,
                              int* passes_out = nullptr)
    {
        using Sys = DeckTrees::Sys;
        const NetworkSolve::Parameters<double> params{1e-2, 80};
        auto fb = dt.build(step, ipr);
        auto outer = dt.build(step, ipr);
        auto reduced = dt.build(step, ipr);
        const auto guess = dt.guess(guess_bar);
        const auto rf = NetworkSolve::solve(fb, guess, params, NetworkSolve::FullStep{});
        const auto ro = NetworkSolve::solveWithTree(outer, guess, params, NetworkSolve::FullStep{});
        const auto rr = NetworkSolve::solveReduced(reduced, guess, params, true);
        std::string sets;
        for (const auto& s : ro.sets) { sets += " " + s; }
        BOOST_TEST_MESSAGE(fmt::format("{}: outer {} passes / {} it{}{} [{}]; reduced {} it, {} stalls;"
                                       " fb {}",
                                       tag, ro.passes, ro.inner_iterations,
                                       ro.consistent ? "" : " NOT CONSISTENT", ro.cycled ? " CYCLED" : "",
                                       sets, rr.iterations, rr.stalls,
                                       rf.converged ? std::to_string(rf.iterations) + " it" : "FAILED"));
        std::string line = "  outer" + rateList(outer, ro.result.well_rate) + "  nodes";
        for (std::size_t n = 1; n < dt.node_order.size(); ++n) {
            line += fmt::format(" {}={:.1f}", dt.node_order[n], convert::to(ro.result.node_pressure[n], bars));
        }
        BOOST_TEST_MESSAGE(line);
        if (rf.converged) {
            std::string f = "  fb   " + rateList(fb, rf.well_rate) + "  nodes";
            for (std::size_t n = 1; n < dt.node_order.size(); ++n) {
                f += fmt::format(" {}={:.1f}", dt.node_order[n], convert::to(rf.node_pressure[n], bars));
            }
            f += "  bhp";
            for (int w = 0; w < fb.numWells(); ++w) {
                f += fmt::format(" {:.1f}", convert::to(rf.well_bhp[w], bars));
            }
            f += "  outer bhp";
            for (int w = 0; w < outer.numWells(); ++w) {
                f += fmt::format(" {:.1f}", convert::to(ro.result.well_bhp[w], bars));
            }
            BOOST_TEST_MESSAGE(f);
        }
        BOOST_CHECK(ro.result.converged);
        BOOST_CHECK(ro.consistent);
        BOOST_CHECK(rr.converged);
        BOOST_CHECK(rf.converged);
        BOOST_CHECK_EQUAL(ro.result.off_axis, 0);
        BOOST_CHECK_EQUAL(rr.off_axis, 0);
        if (passes_out) { *passes_out = ro.passes; }
        if (!ro.result.converged) { return; }
        // Fischer-Burmeister is only warned on: MODEL5's well table has a
        // loading hump, and the continuous well row converges onto its
        // unstable branch -- bhp 199 bar whatever the inflow, three of four
        // configurations -- where the crossing routine takes the stable one.
        for (int w = 0; w < outer.numWells(); ++w) {
            if (rr.converged) { BOOST_CHECK_CLOSE(ro.result.well_rate[w], rr.well_rate[w], 0.5); }
            if (rf.converged) { BOOST_WARN_CLOSE(ro.result.well_rate[w], rf.well_rate[w], 0.5); }
        }
        const auto stein = steinAllocation(outer, ro.result.state, *dt.schedule, step);
        BOOST_CHECK(!stein.empty());
        if (stein.empty()) { return; }
        std::string both;
        double worst = 0.0;
        for (int w = 0; w < outer.numWells(); ++w) {
            const double mine = ro.result.well_rate[w] * 86400.0;
            both += fmt::format(" {}={:.1f}/{:.1f}", outer.wells()[w].name, mine, stein[w]);
            worst = std::max(worst, std::abs(mine - stein[w]) / std::max(std::abs(stein[w]), 1.0));
        }
        BOOST_TEST_MESSAGE("  newton/stein:" << both);
        BOOST_CHECK_LT(worst, 0.005);
    }
}



// ---------------------------------------------------------------------------
// Generated two-tree cases.
//
// Random instances at three sizes, as deck text, so the reader above and
// Stein's oracle apply unchanged. Groups form one tree (GRUPTREE); a subset
// of them are nodes and form another (BRANPROP) with its own parents, so the
// two structures share names and wells and nothing else. Wells sit in random
// groups; targets on about a third of the groups, mostly oil, some liquid,
// set between what the subtree could make and well above it.
// ---------------------------------------------------------------------------
namespace {
    struct GenSpec {
        int wells, nodes, groups; unsigned seed; bool size_tables = true;
        int group_depth = 4, net_depth = 3;      // deepest allowed tree levels
        double stiff = 0.0;                      // fraction of wells with a quarter of the inflow: thp-capped below their share
        double target_fraction = 0.35;           // groups given a target
    };

    std::string vfpprodTable(const int number, const std::vector<double>& flo,
                             const std::vector<double>& thp,
                             const std::function<double(double, double)>& bhp)
    {
        std::string t = fmt::format("VFPPROD\n {} 2000 'LIQ' 'WCT' 'GOR' 'THP' 'GRAT' 'METRIC' 'BHP' /\n", number);
        for (const double f : flo) { t += fmt::format(" {}", f); }
        t += " /\n";
        for (const double p : thp) { t += fmt::format(" {}", p); }
        t += " /\n 0.0 0.5 /\n 50 500 /\n 0 /\n";
        for (std::size_t i = 0; i < thp.size(); ++i) {
            for (int wfr = 1; wfr <= 2; ++wfr) {
                for (int gfr = 1; gfr <= 2; ++gfr) {
                    t += fmt::format(" {} {} {} 1", i + 1, wfr, gfr);
                    for (const double f : flo) { t += fmt::format(" {:.3f}", bhp(thp[i], f)); }
                    t += " /\n";
                }
            }
        }
        return t;
    }

    /// Everything before SCHEDULE: a grid with a cell per well, and the
    /// dimensions the tree and the network need.
    std::string deckPreamble(const int W, const int G, const int N)
    {
        const int nx = 20, ny = (W + nx - 1) / nx, cells = nx * ny * 2;
        std::string d;
        d += "RUNSPEC\nTITLE\nGENERATED_TWO_TREES\nDIMENS\n";
        d += fmt::format(" {} {} 2 /\nOIL\nWATER\nGAS\nDISGAS\nMETRIC\nSTART\n 1 'JAN' 2020 /\n", nx, ny);
        d += fmt::format("WELLDIMS\n {} 2 {} {} /\n", W, G + 2, W);
        d += "EQLDIMS\n 1 1* 25 1* 1 /\nTABDIMS\n 1 1 50 60 1 60 1 1 /\nVFPPDIMS\n 10 8 2 2 1 2 /\n";
        d += fmt::format("NETWORK\n {} {} /\nUNIFOUT\n", N + 1, N + 1);
        d += "GRID\nINIT\nDXV\n" + fmt::format("{}*1000 /\nDYV\n{}*1000 /\nDZV\n50 50 /\nTOPS\n{}*7000 /\n", nx, ny, nx * ny);
        d += fmt::format("PORO\n{0}*0.2 /\nPERMX\n{0}*200 /\nPERMY\n{0}*200 /\nPERMZ\n{0}*20 /\n", cells);
        d += "PROPS\nSWOF\n 0.2 0 1 0\n 0.3 0.07 0.8 0\n 1.0 1 0 0 /\nSGOF\n 0 0 1 0\n 0.05 0 0.8 0\n 0.79 1 0 0 /\n";
        d += "DENSITY\n 800 1000 1 /\nPVTW\n 1 1.0 4.0E-5 0.5 0.0 /\nPVDG\n 1 1.0 0.01\n 100 0.1 0.015\n 300 0.033 0.02 /\n";
        d += "PVTO\n 1 50 1.2 1.0\n 150 1.15 1.1\n 300 1.10 1.2 /\n 10 150 1.25 0.9\n 250 1.20 1.0\n 350 1.15 1.1 /\n/\n";
        d += fmt::format("REGIONS\nEQLNUM\n{}*1 /\nSOLUTION\nEQUIL\n 7000 270 7050 0 7000 0 1* 0 0 /\nSUMMARY\nFOPR\nSCHEDULE\n", cells);
        return d;
    }

    std::string generateDeck(const GenSpec& spec, std::string* summary)
    {
        std::mt19937 rng(spec.seed);
        auto uniform = [&](const double a, const double b) {
            return std::uniform_real_distribution<double>(a, b)(rng);
        };
        auto pick = [&](const int n) { return std::uniform_int_distribution<int>(0, n - 1)(rng); };

        const int W = spec.wells, G = spec.groups, N = spec.nodes;
        const int nx = 20;
        std::string d = deckPreamble(W, G, N);

        // Tables: a well curve with a loading hump, a branch with a linear drop.
        d += vfpprodTable(1, {10, 200, 500, 1000, 1500, 2000, 3000, 4500, 8000, 20000},
                          {5, 10, 20, 40, 60, 80, 100, 130},
                          [](const double thp, const double q) {
                              const double u = (q - 1500.0) / 1500.0;
                              return thp + 90.0 + 30.0 * u * u;
                          });
        // The branch drop scales with the flow the instance will carry, so a
        // root branch under 200 wells stays inside the table's THP axis
        // rather than in its extrapolation.
        const double flow_scale = spec.size_tables ? std::max(1.0, W / 10.0) : 1.0;
        std::vector<double> branch_flo{10, 500, 1000, 2000, 4000, 8000, 12000, 20000, 30000, 50000};
        for (double& f : branch_flo) { f *= flow_scale; }
        d += vfpprodTable(2, branch_flo, {2, 5, 10, 20, 40, 60, 80, 100},
                          [flow_scale](const double thp, const double q) {
                              return thp + 1.0 + 4.0 * q / (5000.0 * flow_scale);
                          });

        // The group tree: PLAT under FIELD, every other group under PLAT or a
        // group before it, depth at most 4.
        std::vector<std::string> gname{"PLAT"};
        std::vector<int> gparent{-1}, gdepth{1};
        for (int g = 1; g < G; ++g) {
            gname.push_back(fmt::format("G{}", g));
            int par;
            do { par = pick(g); } while (gdepth[par] >= spec.group_depth);
            gparent.push_back(par);
            gdepth.push_back(gdepth[par] + 1);
        }
        d += "GRUPTREE\n 'PLAT' 'FIELD' /\n";
        for (int g = 1; g < G; ++g) { d += fmt::format(" '{}' '{}' /\n", gname[g], gname[gparent[g]]); }
        d += "/\n";

        // The flow network: PLAT is the terminal; N-1 other groups are nodes,
        // each hung under an earlier node -- not its group parent.
        std::vector<int> node_of_group(G, -1);
        std::vector<int> nodes{0};
        node_of_group[0] = 0;
        {
            std::vector<int> candidates;
            for (int g = 1; g < G; ++g) { candidates.push_back(g); }
            std::shuffle(candidates.begin(), candidates.end(), rng);
            for (int i = 0; i < std::min(N - 1, G - 1); ++i) {
                node_of_group[candidates[i]] = static_cast<int>(nodes.size());
                nodes.push_back(candidates[i]);
            }
        }
        std::vector<int> ndepth(nodes.size(), 0);
        d += "BRANPROP\n";
        for (std::size_t n = 1; n < nodes.size(); ++n) {
            int up;
            do { up = pick(static_cast<int>(n)); } while (ndepth[up] >= spec.net_depth);
            ndepth[n] = ndepth[up] + 1;
            d += fmt::format(" '{}' '{}' 2 1* /\n", gname[nodes[n]], gname[nodes[up]]);
        }
        d += "/\nNODEPROP\n 'PLAT' 20.0 NO NO 1* /\n";
        for (std::size_t n = 1; n < nodes.size(); ++n) { d += fmt::format(" '{}' 1* NO NO 1* /\n", gname[nodes[n]]); }
        d += "/\n";

        // Wells: a random leaf group -- a group holds wells or groups, not
        // both -- own ORAT limit, bhp limit 100, the well table.
        std::vector<int> leaves;
        for (int g = 0; g < G; ++g) {
            if (std::find(gparent.begin(), gparent.end(), g) == gparent.end()) { leaves.push_back(g); }
        }
        std::vector<int> wgroup(W);
        std::vector<double> worat(W);
        std::string ws = "WELSPECS\n", cd = "COMPDAT\n", wc = "WCONPROD\n";
        for (int w = 0; w < W; ++w) {
            // Every leaf gets a well before any gets a second: an empty group
            // has no guide rate, and the oracle's tree cannot hold one.
            wgroup[w] = (w < static_cast<int>(leaves.size()))
                ? leaves[w] : leaves[pick(static_cast<int>(leaves.size()))];
            worat[w] = std::round(uniform(800.0, 3000.0));
            const int i = w % nx + 1, j = w / nx + 1;
            ws += fmt::format(" 'W{}' '{}' {} {} 7000 'OIL' /\n", w + 1, gname[wgroup[w]], i, j);
            cd += fmt::format(" 'W{}' {} {} 1 2 'OPEN' 1* 1* 0.2 /\n", w + 1, i, j);
            wc += fmt::format(" 'W{}' 'OPEN' 'GRUP' {} 4* 100 1* 1 /\n", w + 1, worat[w]);
        }
        d += ws + "/\n" + cd + "/\n" + wc + "/\n";

        // Targets on about a third of the groups: between 40 % and 90 % of
        // the ORAT limits beneath, or well above them, mostly oil, some liquid.
        std::vector<double> beneath(G, 0.0);
        for (int w = 0; w < W; ++w) {
            for (int g = wgroup[w]; g >= 0; g = gparent[g]) { beneath[g] += worat[w]; }
        }
        std::string gc = "GCONPROD\n";
        int targets = 0, liquid = 0;
        for (int g = 0; g < G; ++g) {
            if (beneath[g] <= 0.0 || uniform(0.0, 1.0) > spec.target_fraction) { continue; }
            const double frac = uniform(0.0, 1.0) < 0.8 ? uniform(0.4, 0.9) : uniform(1.5, 3.0);
            const bool lrat = uniform(0.0, 1.0) < 0.25;
            const double value = std::round(frac * beneath[g] * (lrat ? 1.25 : 1.0));
            gc += lrat ? fmt::format(" '{}' 'LRAT' 1* 1* 1* {} 'RATE' /\n", gname[g], value)
                       : fmt::format(" '{}' 'ORAT' {} 3* 'RATE' /\n", gname[g], value);
            ++targets;
            liquid += lrat ? 1 : 0;
        }
        d += gc + "/\nTSTEP\n1 /\nEND\n";
        if (summary) {
            int ndeep = 0, gdeep = 0;
            for (const int x : ndepth) { ndeep = std::max(ndeep, x); }
            for (const int x : gdepth) { gdeep = std::max(gdeep, x); }
            *summary = fmt::format("{} wells, {} groups ({} with a target, {} of them liquid, depth {}), {} nodes (depth {})",
                                   W, G, targets, liquid, gdeep, nodes.size(), ndeep);
        }
        return d;
    }

    /// The wells a spec means to be weak: every k-th, k from the stiff fraction.
    std::map<std::string, double> weakWells(const GenSpec& spec)
    {
        std::map<std::string, double> j;
        if (spec.stiff <= 0.0) { return j; }
        const int every = std::max(1, static_cast<int>(std::round(1.0 / spec.stiff)));
        for (int w = 0; w < spec.wells; w += every) { j["W" + std::to_string(w + 1)] = 0.25; }
        return j;
    }
}



// A judge that owes nothing to any route: does an answer -- node pressures,
// well rates, the controls it ended on -- satisfy the rows and the
// inequalities that define a solution? Every route's answer on every dumped
// system goes through it, so "best" is decided by the same yardstick.
namespace {
    using NetworkSolve::Verdict;
    using NetworkSolve::verifyAnswer;
}

#endif // OPM_NETWORK_SOLVE_TEST_SUPPORT_HEADER_INCLUDED
