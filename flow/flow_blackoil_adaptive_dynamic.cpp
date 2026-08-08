/*
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

#include <flow/flow_blackoil_adaptive_dynamic.hpp>

#include <opm/material/common/ResetLocale.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicit.hpp>
#include <opm/simulators/flow/FlowMain.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/AdaptiveCpGridVanguard.hpp>
#include <opm/simulators/flow/AdaptiveStateTransfer.hpp>
#include <opm/simulators/flow/python/PyMain.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/Carfin.hpp>
#include <opm/input/eclipse/EclipseState/Grid/LgrCollection.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/models/blackoil/blackoilconvectivemixingmodule.hh>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm::Parameters {

//! \brief Report step at which the simulator world is torn down and rebuilt
//! (grid re-refined from the current mark set, state remapped). -1 = never.
//! Phase-1 exercise of the adaptation seam; the mark set is unchanged for now,
//! so the rebuilt run must reproduce the continuous run.
struct AdaptiveRebuildStep { static constexpr int value = -1; };

//! \brief Refinement spec (same CARFIN-style syntax as --adaptive-lgr) that
//! replaces the mark set at the rebuild step. Empty = keep the initial spec.
struct AdaptiveRebuildLgr { static constexpr auto value = ""; };

} // namespace Opm::Parameters

namespace Opm::Properties {

namespace TTag {
// Same model configuration as FlowProblemAdaptive (flow_blackoil_adaptive.cpp):
// TPFA blackoil with the AdaptiveCpGridVanguard. A separate TypeTag so the two
// executables stay independent while the dynamic driver evolves.
struct FlowProblemAdaptiveDynamic
{ using InheritsFrom = std::tuple<FlowProblem>; };
}

template<class TypeTag>
struct Linearizer<TypeTag, TTag::FlowProblemAdaptiveDynamic>
{ using type = TpfaLinearizer<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlowProblemAdaptiveDynamic>
{ using type = BlackOilLocalResidualTPFA<TypeTag>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowProblemAdaptiveDynamic>
{ static constexpr bool value = false; };

template<class TypeTag>
struct AvoidElementContext<TypeTag, TTag::FlowProblemAdaptiveDynamic>
{ static constexpr bool value = true; };

} // namespace Opm::Properties

namespace Opm {

//! \brief The adaptive vanguard plus the dynamic driver's parameters
//! (registration must happen in the registration phase, so it lives here).
template <class TypeTag>
class AdaptiveDynamicVanguard : public AdaptiveCpGridVanguard<TypeTag>
{
    using Base = AdaptiveCpGridVanguard<TypeTag>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
public:
    explicit AdaptiveDynamicVanguard(Simulator& simulator)
        : Base(simulator)
    {}

    static void registerParameters()
    {
        Base::registerParameters();
        Parameters::Register<Parameters::AdaptiveRebuildStep>(
            "Report step at which the dynamic driver tears down and rebuilds "
            "the simulator (dynamic-refinement phase 1); -1 disables.");
        Parameters::Register<Parameters::AdaptiveRebuildLgr>(
            "Refinement spec applied at the rebuild step (CARFIN-style, see "
            "--adaptive-lgr); empty keeps the initial spec.");
    }

    //! Mark-set override for the rebuild (process-wide: the driver sets it
    //! before constructing the second simulator world).
    static inline std::string rebuildSpecOverride{};

    std::string adaptiveLgrSpec() const override
    {
        if (rebuildSpecOverride == "none") {   // un-mark everything: coarsen
            return {};
        }
        return rebuildSpecOverride.empty()
            ? Base::adaptiveLgrSpec() : rebuildSpecOverride;
    }

    //! Deck-LGR activation for the world being built (LGRON/LGROFF). Names
    //! absent from the map are active, matching ScheduleState::lgr_active().
    static inline std::map<std::string, bool> lgrActivation{};

    static bool lgrIsActive(const std::string& name)
    {
        const auto it = lgrActivation.find(name);
        return (it == lgrActivation.end()) || it->second;
    }

    //! \brief Add deck LGRs honouring LGRON/LGROFF, then any adaptive spec.
    //!
    //! With every deck LGR active (the deck default) this is exactly the base
    //! path. When an LGROFF is in effect only the active subset of the deck
    //! CARFINs is refined, so a rebuild at the LGROFF step coarsens the
    //! switched-off region and a later LGRON rebuild refines it again.
    void addLgrs()
    {
        const auto& lgrs = this->eclState().getLgrs();
        bool anyOff = false;
        for (std::size_t l = 0; l < lgrs.size(); ++l) {
            if (!lgrIsActive(lgrs.getLgr(l).NAME())) {
                anyOff = true;
                break;
            }
        }
        if (!anyOff) {
            Base::addLgrs();
        }
        else {
            std::vector<std::array<int,3>> cellsPerDim, startIJK, endIJK;
            std::vector<std::string> names;
            for (std::size_t l = 0; l < lgrs.size(); ++l) {
                const auto& c = lgrs.getLgr(l);
                if (!lgrIsActive(c.NAME())) {
                    continue;
                }
                if (c.PARENT_NAME() != "GLOBAL") {
                    throw std::runtime_error("LGRON/LGROFF with nested CARFIN ('"
                                             + c.NAME() + "') is not supported");
                }
                cellsPerDim.push_back({c.NX()/(c.I2() + 1 - c.I1()),
                                       c.NY()/(c.J2() + 1 - c.J1()),
                                       c.NZ()/(c.K2() + 1 - c.K1())});
                startIJK.push_back({c.I1(), c.J1(), c.K1()});
                endIJK.push_back({c.I2() + 1, c.J2() + 1, c.K2() + 1});
                names.push_back(c.NAME());
            }
            OpmLog::info("\nLGRON/LGROFF: refining " + std::to_string(names.size())
                         + " of " + std::to_string(lgrs.size()) + " deck LGR(s)");
            if (!names.empty()) {
                this->grid_->addLgrsUpdateLeafView(cellsPerDim, startIJK, endIJK, names);

                // Same post-refinement bookkeeping as the base deck path.
                this->updateGridView_();
                this->updateCartesianToCompressedMapping_();
                this->updateCellDepths_();
                this->updateCellThickness_();
            }
        }

        // Wells and toggled LGRs: give COMPDAT wells an equivalent trajectory
        // (a one-time synthesis; wells that already have one are skipped) and
        // re-derive every trajectory well against the current leaf -- refined
        // or coarse. This puts a well opened together with an LGRON onto the
        // refined child cells, and reverts a well inside a switched-off LGR
        // to its global-grid connections.
        if (lgrs.size() > 0) {
            const auto& inputGrid = this->eclState().getInputGrid();
            this->schedule().synthesizeWellTrajectories(
                [&inputGrid](std::size_t globalIdx)
                { return inputGrid.getCellCenter(globalIdx); },
                [&inputGrid](std::size_t globalIdx)
                { return inputGrid.getCellDims(globalIdx); });

            this->recomputeWellTrajectoriesInLgr_(/*replayOnCoarse=*/true);
        }
    }
};

} // namespace Opm

namespace Opm::Properties {

template<class TypeTag>
struct Vanguard<TypeTag, TTag::FlowProblemAdaptiveDynamic>
{ using type = AdaptiveDynamicVanguard<TypeTag>; };

} // namespace Opm::Properties

namespace Opm {

//! \brief PyMain with a second-init entry: rebuild the FlowMain/simulator
//! world from the already-parsed model description (no re-parse, no second
//! MPI init). One DynamicMain object lives for the whole process.
template <class TypeTag>
class DynamicMain : public PyMain<TypeTag>
{
public:
    using FlowMainType = FlowMain<TypeTag>;
    using PyMain<TypeTag>::PyMain;

    //! Re-populate the vanguard's (moved-from) static model parameters and
    //! construct a fresh FlowMain. The shared_ptrs come from the snapshot the
    //! driver took before the first vanguard construction; the evolving
    //! schedule-state objects are carried over from the previous simulator.
    std::unique_ptr<FlowMainType>
    rebuildFlowBlackoil(const FlowGenericVanguard::SimulationModelParams& snapshot,
                        const Action::State& actionState,
                        const UDQState& udqState)
    {
        FlowGenericVanguard::modelParams_.setupTime_ = snapshot.setupTime_;
        FlowGenericVanguard::modelParams_.eclState_ = snapshot.eclState_;
        FlowGenericVanguard::modelParams_.eclSchedule_ = snapshot.eclSchedule_;
        FlowGenericVanguard::modelParams_.eclSummaryConfig_ = snapshot.eclSummaryConfig_;
        FlowGenericVanguard::modelParams_.actionState_ =
            std::make_unique<Action::State>(actionState);
        FlowGenericVanguard::modelParams_.udqState_ =
            std::make_unique<UDQState>(udqState);
        // Phase-1 limitation: WTEST history restarts empty (the vanguard's
        // WellTestState was moved into the previous well model).
        FlowGenericVanguard::modelParams_.wtestState_ =
            std::make_unique<WellTestState>();

        return flowMainInit<TypeTag>(this->argc_, this->argv_,
                                     this->outputCout_, this->outputFiles_);
    }
};

// ----------------- Main program: the report-step driver -----------------
//
// Instead of FlowMain::execute() (which owns the whole time loop), drive the
// run one report step at a time through the step API the Python bindings use:
//   initFlowBlackoil() -> executeInitStep() -> executeStep()* -> cleanup.
// PyMain is a plain-C++ Main subclass (no Python dependency) that splits
// initialization from the run loop; we reuse it as the bootstrap.
//
// The step loop is the adaptation seam of the dynamic-refinement plan
// (opm-gridrefined docs/DYNAMIC-REFINEMENT-FLOW-PLAN.md, phase 1 steps
// S2-S5): at --adaptive-rebuild-step=N the simulator world is torn down,
// rebuilt from the already-parsed model description (the vanguard re-refines
// during construction), and the extracted state is remapped by stable cell id
// and injected through the restart machinery. With an unchanged mark set the
// rebuilt run must reproduce the continuous run (plan test S8a).
int flowBlackoilTpfaAdaptiveDynamicMainStandalone(int argc, char** argv)
{
    using TypeTag = Properties::TTag::FlowProblemAdaptiveDynamic;

    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    resetLocale();

    // Main's ctor/dtor own MPI_Init/MPI_Finalize.
    auto mainObject = std::make_unique<DynamicMain<TypeTag>>(argc, argv);

    int exitCode = EXIT_SUCCESS;
    auto flowMain = mainObject->initFlowBlackoil(exitCode);
    if (!flowMain) {
        mainObject.reset();
        return exitCode;
    }

    // The vanguard's static model parameters are populated now and will be
    // moved out during the first vanguard construction (inside
    // executeInitStep). Snapshot the shared_ptrs first so a rebuild can
    // re-populate them without re-parsing the deck.
    FlowGenericVanguard::SimulationModelParams snapshot;
    snapshot.setupTime_ = FlowGenericVanguard::modelParams_.setupTime_;
    snapshot.eclState_ = FlowGenericVanguard::modelParams_.eclState_;
    snapshot.eclSchedule_ = FlowGenericVanguard::modelParams_.eclSchedule_;
    snapshot.eclSummaryConfig_ = FlowGenericVanguard::modelParams_.eclSummaryConfig_;

    // LGRON/LGROFF support: the deck's LGR activation at a report step is the
    // schedule's answer to "which CARFINs exist right now". Track it per step;
    // a change is an adaptation event, handled by the same rebuild seam as
    // --adaptive-rebuild-step.
    std::vector<std::string> deckLgrNames;
    {
        const auto& lgrs = snapshot.eclState_->getLgrs();
        for (std::size_t l = 0; l < lgrs.size(); ++l) {
            deckLgrNames.push_back(lgrs.getLgr(l).NAME());
        }
    }
    const auto lgrActivationAt = [&](int step) {
        std::map<std::string, bool> act;
        const auto& sched = *snapshot.eclSchedule_;
        const auto s = std::min<std::size_t>(step, sched.size() - 1);
        for (const auto& name : deckLgrNames) {
            act.emplace(name, sched[s].lgr_active(name));
        }
        return act;
    };
    auto builtActivation = lgrActivationAt(0);
    AdaptiveDynamicVanguard<TypeTag>::lgrActivation = builtActivation;

    int status = flowMain->executeInitStep();
    if (status == EXIT_SUCCESS) {
        const int rebuildStep = Parameters::Get<Parameters::AdaptiveRebuildStep>();

        // One executeStep() per report step. runStep returns a continue flag
        // (false = schedule EXIT), not an exit status; errors throw.
        bool continueLooping = true;
        while (continueLooping && !flowMain->getSimTimer()->done()) {
            const int step = flowMain->getSimTimer()->currentStepNum();
            const auto stepActivation = lgrActivationAt(step);
            if (step == rebuildStep || stepActivation != builtActivation) {
                // ---- the adaptation event (S2-S5, identity mark set) ----
                auto* sim = flowMain->getSimulatorPtr();

                // S3: extract per-cell state keyed by stable cell id, and
                // carry over the evolving schedule-state objects.
                const auto state = extractAdaptiveState<TypeTag>(*sim);
                const Action::State actionState = sim->vanguard().actionState();
                const UDQState udqState = sim->vanguard().udqState();

                // S0b (minimal): a new mark set for world #2, if given.
                if (step == rebuildStep) {
                    const std::string newSpec =
                        Parameters::Get<Parameters::AdaptiveRebuildLgr>();
                    if (!newSpec.empty()) {
                        AdaptiveDynamicVanguard<TypeTag>::rebuildSpecOverride = newSpec;
                    }
                }
                // LGRON/LGROFF: world #2 refines the deck LGRs active now.
                if (stepActivation != builtActivation) {
                    OpmLog::info("LGRON/LGROFF change at report step "
                                 + std::to_string(step) + "; rebuilding the grid");
                    AdaptiveDynamicVanguard<TypeTag>::lgrActivation = stepActivation;
                    builtActivation = stepActivation;
                }

                // S2: tear down world #1, build world #2 from the parsed
                // model description (vanguard refines during construction).
                flowMain.reset();
                flowMain = mainObject->rebuildFlowBlackoil(snapshot, actionState, udqState);
                status = flowMain->executeInitStep();
                if (status != EXIT_SUCCESS) {
                    break;
                }

                // S4+S5: remap onto the new leaf and inject through the
                // restart machinery; position clock, episode and timer.
                auto* sim2 = flowMain->getSimulatorPtr();
                const auto sol = remapAdaptiveState<TypeTag>(state, *sim2);
                injectAdaptiveState<TypeTag>(*sim2, sol, step);
                flowMain->getSimTimer()->setCurrentStepNum(step);
            }
            continueLooping = (flowMain->executeStep() != 0);
        }
        if (status == EXIT_SUCCESS) {
            status = flowMain->executeStepsCleanup();
        }
    }

    flowMain.reset();
    mainObject.reset();   // destructor calls MPI_Finalize
    return status;
}

} // namespace Opm
