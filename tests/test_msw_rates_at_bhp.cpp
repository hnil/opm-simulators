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

#include "config.h"

#define BOOST_TEST_MODULE MswRatesAtBhp

#include "SimulatorFixture.hpp"

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/wells/MultisegmentWell.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cmath>
#include <optional>

namespace Opm::Properties::TTag {
    struct TestMswRatesTypeTag {
        using InheritsFrom = std::tuple<TestTypeTag>;
    };
}

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

// A multisegment well asked for its bhp at the thp limit on a copy of the well state
// (as well testing does) must answer from that copy's segment pressures, not the
// well model's.
BOOST_AUTO_TEST_CASE(BhpAtThpReadsTheGivenWellState)
{
    using TypeTag = Opm::Properties::TTag::TestMswRatesTypeTag;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;
    using MSWell = Opm::MultisegmentWell<TypeTag>;

    auto simulator = Opm::initSimulator<TypeTag>("MSW_THP.DATA", "test_msw_rates_at_bhp", /*threads=*/1);
    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(86400.0);
    simulator->problem().resetIterationForNewTimestep();

    WellModel& well_model = simulator->problem().wellModel();
    well_model.beginReportStep(0);
    well_model.beginTimeStep();
    auto logger_guard = well_model.groupStateHelper().pushLogger();
    auto& deferred_logger = well_model.groupStateHelper().deferredLogger();
    well_model.calculateExplicitQuantities();
    well_model.prepareTimeStep(deferred_logger);

    const auto& well = well_model.getWell("B-1H");
    BOOST_REQUIRE(dynamic_cast<const MSWell*>(&well) != nullptr);
    const auto& summary_state = simulator->vanguard().summaryState();

    auto& global = well_model.wellState();
    const auto original = global;
    auto bhpAtThp = [&](const auto& helper) {
        return well.computeBhpAtThpLimitProdWithAlq(*simulator, helper, summary_state,
                                                    /*alq_value=*/0.0, /*iterate_if_no_solution=*/false);
    };

    // Same top pressure, a different profile below it.
    auto copy = original;
    auto& pressure = copy.well("B-1H").segments.pressure;
    BOOST_REQUIRE(pressure.size() > 2);
    for (std::size_t seg = 1; seg < pressure.size(); ++seg) {
        pressure[seg] += 2.0e5 * seg;
    }

    const auto on_original = bhpAtThp(well_model.groupStateHelper());

    auto helper_copy = well_model.groupStateHelper();
    std::optional<double> on_copy;
    {
        auto guard = helper_copy.pushWellState(copy);
        on_copy = bhpAtThp(helper_copy);
    }

    global = copy;
    const auto copy_as_global = bhpAtThp(well_model.groupStateHelper());
    global = original;

    BOOST_REQUIRE(on_original.has_value());
    BOOST_REQUIRE(on_copy.has_value());
    BOOST_REQUIRE(copy_as_global.has_value());
    // The profile has to matter, or the check below proves nothing.
    BOOST_CHECK_GT(std::abs(*copy_as_global - *on_original), 1.0e3);
    BOOST_CHECK_CLOSE(*on_copy, *copy_as_global, 1.0e-10);
}
