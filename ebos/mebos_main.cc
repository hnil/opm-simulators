// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief The main file of mebos, an multiplexed-version of ebos, the general-purpose
 *        black-oil simulator for ECL decks for research purposes.
 *
 * Just like 'flow', it does not require to select the simulator binary to run a deck
 * that uses certain options like twophase, solvent, polymer or thermal in advance.
 */
#include "config.h"

#include "ebos_blackoil.hh"
#include "ebos_oilwater.hh"
#include "ebos_oilwater_polymer.hh"
#include "ebos_gasoil.hh"
#include "ebos_gaswater.hh"
#include "ebos_energy.hh"
#include "ebos_solvent.hh"
#include "ebos_polymer.hh"
#include "ebos_foam.hh"

#include <opm/models/utils/propertysystem.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Parser/ErrorGuard.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <string>
#include <memory>

int main(int argc, char **argv)
{
    Dune::Timer externalSetupTimer;
    externalSetupTimer.start();

    if (!Opm::ebosBlackOilDeckFileNameIsSet(argc, argv))
        // no deck was specified, e.g., --help. use the black oil variant to figure out
        // what exactly should be done
        return Opm::ebosBlackOilMain(argc, argv);

    std::string deckFileName =
        Opm::ebosBlackOilGetDeckFileName(argc, argv);

    std::unique_ptr<Opm::ParseContext> parseContext
        = Opm::ebosBlackOilCreateParseContext(argc, argv);
    auto errorGuard = std::make_unique<Opm::ErrorGuard>();

    // deal with parallel runs
    int myRank = Dune::MPIHelper::instance(argc, argv).rank();

    Opm::Parser parser;
    // parse the deck file
    if (myRank == 0)
        std::cout << "Parsing deck file \"" << deckFileName << "\"" << std::endl;
    auto deck = std::make_unique<Opm::Deck>(parser.parseFile(deckFileName, *parseContext, *errorGuard));

    // TODO: check which variant ought to be used
    bool waterActive = deck->hasKeyword("WATER");
    bool gasActive = deck->hasKeyword("GAS");
    bool oilActive = deck->hasKeyword("OIL");
    bool solventActive = deck->hasKeyword("SOLVENT");
    bool polymerActive = deck->hasKeyword("POLYMER");
    bool foamActive = deck->hasKeyword("FOAM");
    bool thermalActive = deck->hasKeyword("THERMAL") || deck->hasKeyword("TEMP");

    std::stringstream notSupportedErrorStream;
    notSupportedErrorStream << "deck not supported by mebos, you might want to use a specialized binary. Active options:\n"
                            << "   water: " << waterActive << "\n"
                            << "   gas: " << gasActive << "\n"
                            << "   oil: " << oilActive << "\n"
                            << "   solvent: " << solventActive << "\n"
                            << "   polymer: " << polymerActive << "\n"
                            << "   foam: " << foamActive << "\n"
                            << "   thermal/temperature: " << thermalActive << "\n";

    int numBlackOilPhases = (waterActive?1:0) + (gasActive?1:0) + (oilActive?1:0);
    if (numBlackOilPhases == 0) {
        notSupportedErrorStream << "\n"
                                << "no black-oil phase (water, gas or oil) specified.\n";
        std::cerr << notSupportedErrorStream.str() << std::endl;
        std::abort();
    }
    else if (numBlackOilPhases == 1) {
        notSupportedErrorStream << "\n"
                                << "single-phase simulations are unsupported\n";
        std::cerr << notSupportedErrorStream.str() << std::endl;
        std::abort();
    }
    else if (numBlackOilPhases == 2) {
        if (solventActive) {
            notSupportedErrorStream << "\n"
                                    << "combining twophase and solvent is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (polymerActive && oilActive && waterActive) {
            if (myRank == 0)
                std::cout << "Using oil-water-polymer mode" << std::endl;
            Opm::ebosOilWaterPolymerSetDeck(std::move(deck),
                                            std::move(parseContext),
                                            std::move(errorGuard),
                                            externalSetupTimer.elapsed());
            return Opm::ebosOilWaterPolymerMain(argc, argv);
        }

        if (polymerActive) {
            notSupportedErrorStream << "\n"
                                    << "combining twophase and polymer is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (foamActive) {
            notSupportedErrorStream << "\n"
                                    << "combining twophase and foam is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (thermalActive) {
            notSupportedErrorStream << "\n"
                                    << "combining twophase and energy conservation is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (oilActive && waterActive) {
            if (myRank == 0)
                std::cout << "Using oil-water mode" << std::endl;
            Opm::ebosOilWaterSetDeck(std::move(deck),
                                     std::move(parseContext),
                                     std::move(errorGuard),
                                     externalSetupTimer.elapsed());
            return Opm::ebosOilWaterMain(argc, argv);
        }
        else if (oilActive && gasActive) {
            // run ebos_gasoil
            if (myRank == 0)
                std::cout << "Using gas-oil mode" << std::endl;
            Opm::ebosGasOilSetDeck(std::move(deck),
                                   std::move(parseContext),
                                   std::move(errorGuard),
                                   externalSetupTimer.elapsed());
            return Opm::ebosGasOilMain(argc, argv);
        }
        else if (waterActive && gasActive) {
            // run ebos_gaswater
            if (myRank == 0)
                std::cout << "Using gas-water mode" << std::endl;
            Opm::ebosGasWaterSetDeck(std::move(deck),
                                     std::move(parseContext),
                                     std::move(errorGuard),
                                     externalSetupTimer.elapsed());
            return Opm::ebosGasWaterMain(argc, argv);
        }
    }
    else if (foamActive) {
        if (solventActive) {
            notSupportedErrorStream << "\n"
                                    << "combining foam and solvent is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (polymerActive) {
            notSupportedErrorStream << "\n"
                                    << "combining foam and polymer is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (thermalActive) {
            notSupportedErrorStream << "\n"
                                    << "combining foam and and energy conservation is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        // run ebos_foam
        if (myRank == 0)
            std::cout << "Using foam mode" << std::endl;
        Opm::ebosFoamSetDeck(std::move(deck),
                             std::move(parseContext),
                             std::move(errorGuard),
                             externalSetupTimer.elapsed());
        return Opm::ebosFoamMain(argc, argv);
    }
    else if (polymerActive) {
        if (solventActive) {
            notSupportedErrorStream << "\n"
                                    << "combining polymer and solvent is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (foamActive) {
            notSupportedErrorStream << "\n"
                                    << "combining polymer and foam is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (thermalActive) {
            notSupportedErrorStream << "\n"
                                    << "combining polymer and and energy conservation is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        // run ebos_polymer
        if (myRank == 0)
            std::cout << "Using polymer mode" << std::endl;
        Opm::ebosPolymerSetDeck(std::move(deck),
                                std::move(parseContext),
                                std::move(errorGuard),
                                externalSetupTimer.elapsed());
        return Opm::ebosPolymerMain(argc, argv);
    }
    else if (solventActive) {
        if (polymerActive) {
            notSupportedErrorStream << "\n"
                                    << "combining solvent and polymer is not supported\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (foamActive) {
            notSupportedErrorStream << "\n"
                                    << "combining solvent and foam is not supported\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (thermalActive) {
            notSupportedErrorStream << "\n"
                                    << "combining solvent and and energy conservation is not supported\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        // run ebos_solvent
        if (myRank == 0)
            std::cout << "Using solvent mode" << std::endl;
        Opm::ebosSolventSetDeck(std::move(deck),
                                std::move(parseContext),
                                std::move(errorGuard),
                                externalSetupTimer.elapsed());
        return Opm::ebosSolventMain(argc, argv);
    }
    else if (thermalActive) {
        if (solventActive) {
            notSupportedErrorStream << "\n"
                                    << "combining thermal and solvent is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (polymerActive) {
            notSupportedErrorStream << "\n"
                                    << "combining thermal and polymer is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (foamActive) {
            notSupportedErrorStream << "\n"
                                    << "combining thermal and foam is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        // run ebos_thermal
        if (myRank == 0)
            std::cout << "Using thermal mode" << std::endl;
        Opm::ebosEnergySetDeck(std::move(deck),
                               std::move(parseContext),
                               std::move(errorGuard),
                               externalSetupTimer.elapsed());
        return Opm::ebosEnergyMain(argc, argv);
    }
    else {
        if (myRank == 0)
            std::cout << "Using blackoil mode" << std::endl;
        Opm::ebosBlackOilSetDeck(std::move(deck),
                                 std::move(parseContext),
                                 std::move(errorGuard),
                                 externalSetupTimer.elapsed());
        return Opm::ebosBlackOilMain(argc, argv);
    }

    if (myRank == 0)
        // this is supposed to be unreachable. this should not happen!
        std::cerr << "Oops: something went wrong when deciding which simulator ought to be used" << std::endl;
    std::abort();

    return 0;
}