/*
  Copyright 2026 Equinor ASA

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

#define BOOST_TEST_MODULE OPM_test_FluxActivation

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/simulators/flow/flux/FluxActivation.hpp>

#include <stdexcept>
#include <string>

namespace {

Opm::Deck makeDeck(const std::string& gridKeywords)
{
    const auto data = std::string{R"(
RUNSPEC
OIL
METRIC
DIMENS
 3 1 1 /
TABDIMS
 1 1 40 20 1 20 /
GRID
)"}
        + gridKeywords
        + R"(
DX
 3*1 /
DY
 3*1 /
DZ
 3*1 /
TOPS
 3*0 /
PROPS
ROCK
 14.7 3E-6 /
PVDO
 14.7 1.062 1.04 /
DENSITY
 53.66 64.49 0.0533 /
SOLUTION
PRESSURE
 3*300 /
SCHEDULE
TSTEP
 1 /
END
 )";

    return Opm::Parser{}.parseString(data);
}

} // namespace

BOOST_AUTO_TEST_CASE(DoesNothingWithoutUseFlux)
{
    auto deck = makeDeck("FLUXNUM\n 1 0 1 /\n");
    Opm::EclipseState eclipseState(deck);

    const auto originalActnum = eclipseState.fieldProps().actnum();
    BOOST_CHECK(!Opm::applyUseFluxActnum(eclipseState));
    BOOST_CHECK(eclipseState.fieldProps().actnum() == originalActnum);
}

BOOST_AUTO_TEST_CASE(ActivatesOnlyTheSelectedFluxnumRegion)
{
    auto deck = makeDeck("FLUXNUM\n 1 0 1 /\nUSEFLUX\n 1* /\n");
    Opm::EclipseState eclipseState(deck);

    BOOST_CHECK(Opm::applyUseFluxActnum(eclipseState));

    const auto actnum = eclipseState.fieldProps().actnum();
    BOOST_REQUIRE_EQUAL(actnum.size(), 3U);
    BOOST_CHECK_EQUAL(actnum[0], 1);
    BOOST_CHECK_EQUAL(actnum[1], 0);
    BOOST_CHECK_EQUAL(actnum[2], 1);

    const auto& gridActnum = eclipseState.getInputGrid().getACTNUM();
    BOOST_REQUIRE_EQUAL(gridActnum.size(), 3U);
    BOOST_CHECK_EQUAL(gridActnum[0], 1);
    BOOST_CHECK_EQUAL(gridActnum[1], 0);
    BOOST_CHECK_EQUAL(gridActnum[2], 1);
}

BOOST_AUTO_TEST_CASE(DoesNotReactivatePreviouslyInactiveCells)
{
    auto deck = makeDeck("ACTNUM\n 1 0 1 /\nFLUXNUM\n 1 1 1 /\nFLUXREG\n 1 /\nUSEFLUX\n 1* /\n");
    Opm::EclipseState eclipseState(deck);

    BOOST_CHECK(Opm::applyUseFluxActnum(eclipseState));

    const auto actnum = eclipseState.fieldProps().actnum();
    BOOST_REQUIRE_EQUAL(actnum.size(), 3U);
    BOOST_CHECK_EQUAL(actnum[0], 1);
    BOOST_CHECK_EQUAL(actnum[1], 0);
    BOOST_CHECK_EQUAL(actnum[2], 1);
}

BOOST_AUTO_TEST_CASE(RejectsUseFluxWithoutFluxnum)
{
    auto deck = makeDeck("FLUXREG\n 1 /\nUSEFLUX\n 1* /\n");
    Opm::EclipseState eclipseState(deck);

    BOOST_CHECK_THROW(Opm::applyUseFluxActnum(eclipseState), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(ActivatesSelectedFluxregIds)
{
    auto deck = makeDeck("FLUXNUM\n 1 2 3 /\nFLUXREG\n 1 3 /\nUSEFLUX\n 1* /\n");
    Opm::EclipseState eclipseState(deck);

    BOOST_CHECK(Opm::applyUseFluxActnum(eclipseState));

    const auto actnum = eclipseState.fieldProps().actnum();
    BOOST_REQUIRE_EQUAL(actnum.size(), 3U);
    BOOST_CHECK_EQUAL(actnum[0], 1);
    BOOST_CHECK_EQUAL(actnum[1], 0);
    BOOST_CHECK_EQUAL(actnum[2], 1);
}

BOOST_AUTO_TEST_CASE(RejectsAmbiguousFluxnumSelection)
{
    auto deck = makeDeck("FLUXNUM\n 1 2 0 /\nUSEFLUX\n 1* /\n");
    Opm::EclipseState eclipseState(deck);

    BOOST_CHECK_THROW(Opm::applyUseFluxActnum(eclipseState), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(RejectsFluxregSelectingUnknownFluxnumRegion)
{
    auto deck = makeDeck("FLUXNUM\n 1 2 0 /\nFLUXREG\n 3 /\nUSEFLUX\n 1* /\n");
    Opm::EclipseState eclipseState(deck);

    BOOST_CHECK_THROW(Opm::applyUseFluxActnum(eclipseState), std::invalid_argument);
}
