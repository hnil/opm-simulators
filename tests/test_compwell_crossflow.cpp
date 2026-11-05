/*
  Copyright 2026, SINTEF Digital

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
#include "config.h"

#define BOOST_TEST_MODULE CompWellCrossflow
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>

#include <flowexperimental/comp/flow_comp.hpp>

#include "RunDeckHelper.hpp"

#include <opm/output/eclipse/WStat.hpp>

#include <cstddef>
#include <string>

namespace {

// Two layers without flow between them, at 100 and 60 bar, and a producer
// at 70 bar connected to both. With crossflow, the producer injects into the
// lower-pressure layer through its connection there.
// drain_upper_layer adds a second producer on the 100 bar layer with a 65 bar limit.
std::string deck(const std::string& crossflow,
                 const std::string& bhp_limit = "70",
                 const bool drain_upper_layer = false)
{
    const std::string drain_welspecs = drain_upper_layer ? "DRAIN FIELD 1 1 1* GAS /\n" : "";
    const std::string drain_compdat = drain_upper_layer ? "DRAIN 1 1 1 1 OPEN 2* 0.0151 /\n" : "";
    const std::string drain_wconprod = drain_upper_layer ? "DRAIN OPEN BHP 2* 100000 2* 65 /\n" : "";
    const std::string steps = drain_upper_layer ? "20*0.1" : "5*0.1";
    return R"(RUNSPEC
METRIC
TABDIMS
/
OIL
GAS
DIMENS
1 1 2 /
COMPS
3 /
FULLIMP
WELLDIMS
2 2 1 2 /
START
1 'JAN' 2016 /
UNIFOUT

GRID
DX
2*1000 /
DY
2*1 /
DZ
2*1 /
TOPS
0 /
PERMX
2*100 /
PERMY
2*100 /
PERMZ
2*0 /
PORO
2*0.1 /

PROPS
CNAMES
CO2
METHANE
DECANE
/
ROCK
68.9476 0 /
EOS
PR /
BIC
0
0
0
/
ACF
0.22394 0.01142 0.4884 /
PCRIT
73.773 45.992 21.03 /
TCRIT
304.128 190.564 617.7 /
MW
44.00 16.04 142.28 /
VCRIT
0.09412 0.09863 0.60980 /
SGOF
0.0 0.00 1.00 0.0
0.2 0.04 0.64 0.0
0.4 0.16 0.36 0.0
0.6 0.36 0.16 0.0
0.8 0.64 0.04 0.0
1.0 1.00 0.00 0.0
/
STCOND
15.0 1.0 /

SOLUTION
PRESSURE
100 60 /
SGAS
2*1 /
TEMPI
2*150 /
ZMF
2*0.1
2*0.3
2*0.6
/

SUMMARY
WGPR
/
WOPR
/
WSTAT
/
BPR
1 1 1 /
1 1 2 /
/

SCHEDULE
WELSPECS
PROD FIELD 1 1 1* GAS 3* )" + crossflow + R"( /
)" + drain_welspecs + R"(/
COMPDAT
PROD 1 1 1 2 OPEN 2* 0.0151 /
)" + drain_compdat + R"(/
WCONPROD
PROD OPEN BHP 2* 100000 2* )" + bhp_limit + R"( /
)" + drain_wconprod + R"(/
TSTEP
)" + steps + R"( /
END
)";
}

// The deck has no water phase.
using TypeTag = Opm::Properties::TTag::FlowCompProblem<3, false>;

void checkStoppedAtEnd(const Opm::EclIO::ESmry& smry)
{
    BOOST_CHECK_EQUAL(smry.get("WSTAT:PROD").back(),
                      static_cast<float>(Opm::WStat::numeric::STOP));
    BOOST_CHECK_SMALL(smry.get("WGPR:PROD").back(), 1.0e-6f);
    BOOST_CHECK_SMALL(smry.get("WOPR:PROD").back(), 1.0e-6f);
    BOOST_CHECK_SMALL(smry.get("BPR:1,1,2").back() - 60.0f, 1.0e-2f);
}

} // namespace

using RunDeckFixture = Opm::RunDeckFixture;
BOOST_GLOBAL_FIXTURE(RunDeckFixture);

// Without crossflow, the producer's connection in the lower-pressure layer is
// closed, and that layer keeps its pressure, while the producer still produces
// from the other layer. With crossflow, the layer takes up fluid.
BOOST_AUTO_TEST_CASE(ProducerConnectionDoesNotInjectWithoutCrossflow)
{
    // Give each run a separate output directory.
    const Opm::TemporaryDirectory dir("compwell_crossflow");

    const auto without = Opm::runDeck<TypeTag>(dir.path(), "WITHOUT_CROSSFLOW", deck("NO"));
    const auto& time = without.get("TIME");
    BOOST_CHECK_GT(time.size(), 0u);
    for (std::size_t step = 0; step < time.size(); ++step) {
        BOOST_TEST_CONTEXT("without crossflow at day " << time[step]) {
            BOOST_CHECK_SMALL(without.get("BPR:1,1,2")[step] - 60.0f, 1.0e-3f);
            BOOST_CHECK_LT(without.get("BPR:1,1,1")[step], 100.0f);
            BOOST_CHECK_GT(without.get("WGPR:PROD")[step], 1.0f);
        }
    }

    const auto with = Opm::runDeck<TypeTag>(dir.path(), "WITH_CROSSFLOW", deck("YES"));
    BOOST_CHECK_GT(with.get("BPR:1,1,2").back(), 61.0f);
}

// Both layers are below the producer's BHP limit, so it stops: no surface rate
// and, without crossflow, no flow into either layer.
BOOST_AUTO_TEST_CASE(ProducerStopsWithoutCrossflowWhenAboveAllLayers)
{
    const Opm::TemporaryDirectory dir("compwell_crossflow_stop");

    const auto smry = Opm::runDeck<TypeTag>(dir.path(), "ABOVE_ALL_LAYERS", deck("NO", "110"));
    checkStoppedAtEnd(smry);
    BOOST_CHECK_SMALL(smry.get("BPR:1,1,1").back() - 100.0f, 1.0e-2f);
}

// PROD stops once the 100 bar layer, which DRAIN depletes, falls below PROD's 70 bar limit.
BOOST_AUTO_TEST_CASE(ProducerStopsWithoutCrossflowWhenItsLayerDepletes)
{
    const Opm::TemporaryDirectory dir("compwell_crossflow_deplete");

    const auto smry = Opm::runDeck<TypeTag>(dir.path(), "DEPLETING_LAYER", deck("NO", "70", true));
    BOOST_CHECK_GT(smry.get("WGPR:PROD").front(), 1.0f);
    BOOST_CHECK_LT(smry.get("BPR:1,1,1").back(), 70.0f);
    checkStoppedAtEnd(smry);
}
