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
 * \copydoc Opm::EclWriter
 */
#ifndef OPM_ECL_WRITER_HPP
#define OPM_ECL_WRITER_HPP

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/partitionset.hh>

#include <opm/common/TimingMacros.hpp> // OPM_TIMEBLOCK
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/TimeService.hpp>
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionX.hpp>
#include <opm/input/eclipse/Schedule/RPTConfig.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>

#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/io/eclipse/SummaryNode.hpp>

#include <opm/output/eclipse/Inplace.hpp>
#include <opm/output/eclipse/RegionVariableCollection.hpp>
#include <opm/output/eclipse/RestartValue.hpp>

#include <opm/models/blackoil/blackoilproperties.hh> // Properties::EnableMech, EnableSolvent
#include <opm/models/common/multiphasebaseproperties.hh> // Properties::FluidSystem

#include <opm/simulators/flow/CollectDataOnIORank.hpp>
#include <opm/simulators/flow/countGlobalCells.hpp>
#include <opm/simulators/flow/EclGenericWriter.hpp>
#include <opm/simulators/flow/FlowBaseVanguard.hpp>
#include <opm/simulators/flow/FlowProblemParameters.hpp>
#include <opm/simulators/flow/flux/FluxDumper.hpp>
#include <opm/simulators/flow/flux/FluxParentWells.hpp>
#include <opm/simulators/flow/flux/ParentSummary.hpp>
#include <opm/simulators/flow/flux/FluxSummaryKeys.hpp>
#include <opm/simulators/flow/flux/FluxRegions.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelRestart.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>

#include <opm/simulators/flow/rescoup/ReservoirCouplingEnabled.hpp>
#ifdef RESERVOIR_COUPLING_ENABLED
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#endif

#include <boost/date_time/posix_time/posix_time.hpp>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <cstdio>
#include <unordered_set>
#include <utility>
#include <vector>

namespace Opm::Parameters {

// If available, write the ECL output in a non-blocking manner
struct EnableAsyncEclOutput { static constexpr bool value = true; };

// By default, use single precision for the ECL formated results
struct EclOutputDoublePrecision { static constexpr bool value = false; };

// Write all solutions for visualization, not just the ones for the
// report steps...
struct EnableWriteAllSolutions { static constexpr bool value = false; };

// Write ESMRY file for fast loading of summary data
struct EnableEsmry { static constexpr bool value = true; };

} // namespace Opm::Parameters

namespace Opm::Action {
    class State;
} // namespace Opm::Action

namespace Opm {
    class EclipseIO;
    class UDQState;
} // namespace Opm

namespace Opm {
/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Collects necessary output values and pass it to opm-common's ECL output.
 *
 * Caveats:
 * - For this class to do do anything meaningful, you will have to
 *   have the OPM module opm-common with ECL writing enabled.
 * - The only DUNE grid which is currently supported is Dune::CpGrid
 *   from the OPM module "opm-grid". Using another grid won't
 *   fail at compile time but you will provoke a fatal exception as
 *   soon as you try to write an ECL output file.
 * - This class requires to use the black oil model with the element
 *   centered finite volume discretization.
 */
template <class TypeTag, class OutputModule>
class EclWriter : public EclGenericWriter<GetPropType<TypeTag, Properties::Grid>,
                                          GetPropType<TypeTag, Properties::EquilGrid>,
                                          GetPropType<TypeTag, Properties::GridView>,
                                          GetPropType<TypeTag, Properties::ElementMapper>,
                                          GetPropType<TypeTag, Properties::Scalar>>
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;
    using BaseType = EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>;

    typedef Dune::MultipleCodimMultipleGeomTypeMapper< GridView > VertexMapper;

    static constexpr bool enableEnergy =
        getPropValue<TypeTag, Properties::EnergyModuleType>() == EnergyModules::FullyImplicitThermal ||
        getPropValue<TypeTag, Properties::EnergyModuleType>() == EnergyModules::SequentialImplicitThermal;
    enum { enableMech = getPropValue<TypeTag, Properties::EnableMech>() };
    static constexpr bool enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>();
    enum { enableGeochemistry = getPropValue<TypeTag, Properties::EnableGeochemistry>() };

public:

    using DynamicConns =
        std::vector<std::pair<std::string, std::vector<std::size_t>>>;

    static void registerParameters()
    {
        OutputModule::registerParameters();

        Parameters::Register<Parameters::EnableAsyncEclOutput>
            ("Write the ECL-formated results in a non-blocking way "
             "(i.e., using a separate thread).");
        Parameters::Register<Parameters::EnableEsmry>
            ("Write ESMRY file for fast loading of summary data.");
        Parameters::Register<Parameters::FluxSummaryMinIntervalBetweenSamples<Scalar>>
            ("Minimum time in days between consecutive parent summary samples "
             "written to the .FLUX file by a DUMPFLUX run. Every time step is "
             "sampled unless that would place the sample within this interval "
             "of the previous one; report step boundaries are always sampled. "
             "Smaller values improve the accuracy of the interpolation done by "
             "a USEFLUX run, at the cost of a larger .FLUX file. Use 0 to "
             "sample every time step.");
        Parameters::Register<Parameters::FluxBoundaryMinIntervalBetweenSamples<Scalar>>
            ("Minimum time in days between consecutive sector boundary records "
             "written to the .FLUX file by a DUMPFLUX run. Every time step is "
             "written unless that would place the record within this interval "
             "of the previous one; report step boundaries are always written. "
             "Use 0 to write every time step. Ignored when "
             "--flux-boundary-report-steps-only is set.");
        Parameters::Register<Parameters::FluxBoundaryReportStepsOnly>
            ("Write sector boundary data once per report step rather than per "
             "time step. This restores the behaviour of earlier versions and is "
             "generally less accurate for models whose boundary flow varies "
             "within a report step.");
        Parameters::Register<Parameters::FluxBoundaryInterpolate>
            ("In a USEFLUX run, interpolate the pressure-mode boundary state "
             "(pressure, saturations, Rs, Rv, relative permeability, capillary "
             "pressure) linearly in time between consecutive .FLUX records, "
             "evaluated at the end of each time step. When false, each record's "
             "state is held over the whole interval it covers. The two agree "
             "wherever this run's steps end on the records' times.");
    }

    // The Simulator object should preferably have been const - the
    // only reason that is not the case is due to the SummaryState
    // object owned deep down by the vanguard.
    explicit EclWriter(Simulator& simulator)
        : BaseType(simulator.vanguard().schedule(),
                   simulator.vanguard().eclState(),
                   simulator.vanguard().summaryConfig(),
                   simulator.vanguard().grid(),
                   ((simulator.vanguard().grid().comm().rank() == 0)
                    ? &simulator.vanguard().equilGrid()
                    : nullptr),
                   simulator.vanguard().gridView(),
                   simulator.vanguard().cartesianIndexMapper(),
                   ((simulator.vanguard().grid().comm().rank() == 0)
                    ? &simulator.vanguard().equilCartesianIndexMapper()
                    : nullptr),
                   Parameters::Get<Parameters::EnableAsyncEclOutput>(),
                   Parameters::Get<Parameters::EnableEsmry>())
        , simulator_(simulator)
    {
#if HAVE_MPI
        if (this->simulator_.vanguard().grid().comm().size() > 1) {
            auto smryCfg = (this->simulator_.vanguard().grid().comm().rank() == 0)
                ? this->eclIO_->finalSummaryConfig()
                : SummaryConfig{};

            eclBroadcast(this->simulator_.vanguard().grid().comm(), smryCfg);

            this->outputModule_ = std::make_unique<OutputModule>
                (simulator, smryCfg, this->collectOnIORank_);
        }
        else
#endif
        {
            this->outputModule_ = std::make_unique<OutputModule>
                (simulator, this->eclIO_->finalSummaryConfig(), this->collectOnIORank_);
        }

        this->rank_ = this->simulator_.vanguard().grid().comm().rank();
        this->checkWellsWithinSingleFluxRegion_();
        this->initializeFluxDumpers_();

        this->simulator_.vanguard().eclState().computeFipRegionStatistics();
    }

    ~EclWriter()
    {}

    const EquilGrid& globalGrid() const
    {
        return simulator_.vanguard().equilGrid();
    }

    void recordNewDynamicWellConns(const DynamicConns& newConns)
    {
        if (this->collectOnIORank_.isIORank() && (this->eclIO_ != nullptr)) {
            this->eclIO_->recordNewDynamicWellConns(newConns);
        }
    }

    /*!
     * \brief collect and pass data and pass it to eclIO writer
     */
    void evalSummaryState(bool isSubStep)
    {
        OPM_TIMEBLOCK(evalSummaryState);
        const int reportStepNum = simulator_.episodeIndex() + 1;

        /*
          The summary data is not evaluated for timestep 0, that is
          implemented with a:

             if (time_step == 0)
                 return;

          check somewhere in the summary code. When the summary code was
          split in separate methods Summary::eval() and
          Summary::add_timestep() it was necessary to pull this test out
          here to ensure that the well and group related keywords in the
          restart file, like XWEL and XGRP were "correct" also in the
          initial report step.

          "Correct" in this context means unchanged behavior, might very
          well be more correct to actually remove this if test.
        */

        if (reportStepNum == 0)
            return;

        const Scalar curTime = simulator_.time() + simulator_.timeStepSize();
        const Scalar totalCpuTime =
            simulator_.executionTimer().realTimeElapsed() +
            simulator_.setupTimer().realTimeElapsed() +
            simulator_.vanguard().setupTime();

        auto& regVars = this->outputModule_->regionVariables();

        regVars.prepareValueAccumulation();

        if (const auto conn_opt_ix = regVars
            .variableIndex(this->outputModule_->regVarMapping(), "ConnOPT");
            conn_opt_ix.has_value())
        {
            this->simulator_.problem()
                .wellModel().reportIntervalConnectionOilProduction
                (this->simulator_.timeStepSize(), *conn_opt_ix, regVars);
        }

        auto localWellData                  = simulator_.problem().wellModel().wellData();
        const auto localWBP                 = simulator_.problem().wellModel().wellBlockAveragePressures();
        const auto localGroupAndNetworkData = simulator_.problem().wellModel()
            .groupAndNetworkData(reportStepNum);

        const auto localAquiferData = simulator_.problem().aquiferModel().aquiferData();
        const auto localWellTestState = simulator_.problem().wellModel().wellTestState();
        this->prepareLocalCellData(isSubStep, reportStepNum);

        // A sector run reports FPR over its own cells. FPRFLUX covers the
        // whole of the model the sector came from, using the sums its parent
        // recorded for everything outside.
        {
            std::array<Scalar, 8> hydrocarbonPvWeighted{};
            std::array<Scalar, 8> poreVolumeWeighted{};
            std::array<Scalar, 16> sums{};
            if (this->simulator_.problem()
                    .fluxConverterExternalSums(hydrocarbonPvWeighted, poreVolumeWeighted))
            {
                std::copy(hydrocarbonPvWeighted.begin(), hydrocarbonPvWeighted.end(),
                          sums.begin());
                std::copy(poreVolumeWeighted.begin(), poreVolumeWeighted.end(),
                          sums.begin() + 8);
            }

            this->outputModule_->setExternalRegionSums(sums);
        }

        if (this->outputModule_->needInterfaceFluxes(isSubStep)) {
            this->captureLocalFluxData();
        }

        if (this->collectOnIORank_.isParallel()) {
            OPM_BEGIN_PARALLEL_TRY_CATCH()

            std::map<std::pair<std::string,int>,double> dummy;
            this->collectOnIORank_.collect({},
                                           outputModule_->getBlockData(),
                                           dummy,
                                           localWellData,
                                           localWBP,
                                           localGroupAndNetworkData,
                                           localAquiferData,
                                           localWellTestState,
                                           this->outputModule_->getInterRegFlows(),
                                           {},
                                           {},
                                           this->outputModule_->getLgrBlockData());

            if (this->collectOnIORank_.isIORank()) {
                auto& iregFlows = this->collectOnIORank_.globalInterRegFlows();

                if (! iregFlows.readIsConsistent()) {
                    throw std::runtime_error {
                        "Inconsistent inter-region flow "
                        "region set names in parallel"
                    };
                }

                iregFlows.compress();
            }

            OPM_END_PARALLEL_TRY_CATCH("Collect to I/O rank: ",
                                       this->simulator_.vanguard().grid().comm());
        }


        std::map<std::string, double> miscSummaryData;
        std::map<std::string, std::vector<double>> regionData;
        Inplace inplace;

        {
            OPM_TIMEBLOCK(outputFipLogAndFipresvLog);

            inplace = outputModule_->calc_inplace(miscSummaryData, regionData, simulator_.gridView().comm());

            if (this->collectOnIORank_.isIORank()){
                inplace_ = inplace;
            }
        }

        // Add TCPU
        if (totalCpuTime != 0.0) {
            miscSummaryData["TCPU"] = totalCpuTime;
        }
        if (this->sub_step_report_.total_newton_iterations != 0) {
            miscSummaryData["NEWTON"] = this->sub_step_report_.total_newton_iterations;
        }
        if (this->sub_step_report_.total_linear_iterations != 0) {
            miscSummaryData["MLINEARS"] = this->sub_step_report_.total_linear_iterations;
        }
        if (this->sub_step_report_.total_newton_iterations != 0) {
            miscSummaryData["NLINEARS"] =  static_cast<float>(this->sub_step_report_.total_linear_iterations) / this->sub_step_report_.total_newton_iterations;
        }
        if (this->sub_step_report_.min_linear_iterations != std::numeric_limits<unsigned int>::max()) {
            miscSummaryData["NLINSMIN"] = this->sub_step_report_.min_linear_iterations;
        }
        if (this->sub_step_report_.max_linear_iterations != 0) {
            miscSummaryData["NLINSMAX"] = this->sub_step_report_.max_linear_iterations;
        }
        if (this->simulation_report_.success.total_linear_iterations != 0) {
            miscSummaryData["MSUMLINS"] = this->simulation_report_.success.total_linear_iterations;
        }
        if (this->simulation_report_.success.total_newton_iterations != 0) {
            miscSummaryData["MSUMNEWT"] = this->simulation_report_.success.total_newton_iterations;
        }

        // For reservoir coupling master: collect slave production/injection
        // rates to pass through to Summary::eval() via DynamicSimulatorState.
        const auto rcGroupRates = this->collectReservoirCouplingGroupRates_();

        {
            OPM_TIMEBLOCK(evalSummary);

            // Note: This statement sums one value per registered region
            // variable per region per registered region set across all MPI
            // ranks.
            regVars.commitValues();

            // Wells outside the USEFLUX region take no part in the reduced run,
            // so without this the field and group aggregates would only cover
            // the sector. Filling in their rates from the parent run lets
            // Summary::eval() form F* and G* vectors over the union of sector
            // and parent-only wells through its normal code path.
            this->injectParentOnlyWellData_(this->collectOnIORank_.isParallel()
                                            ? this->collectOnIORank_.globalWellData()
                                            : localWellData,
                                            curTime);

            const auto& blockData = this->collectOnIORank_.isParallel()
                ? this->collectOnIORank_.globalBlockData()
                : this->outputModule_->getBlockData();

            const auto& lgrBlockData = this->collectOnIORank_.isParallel()
                ? this->collectOnIORank_.globalLgrBlockData()
                : this->outputModule_->getLgrBlockData();

            const auto& interRegFlows = this->collectOnIORank_.isParallel()
                ? this->collectOnIORank_.globalInterRegFlows()
                : this->outputModule_->getInterRegFlows();

            this->evalSummary(reportStepNum,
                              curTime,
                              localWellData,
                              localWBP,
                              localGroupAndNetworkData,
                              localAquiferData,
                              blockData,
                              lgrBlockData,
                              miscSummaryData,
                              regionData,
                              this->outputModule_->regVarMapping(),
                              regVars,
                              inplace,
                              this->outputModule_->initialInplace(),
                              interRegFlows,
                              this->summaryState(),
                              this->udqState(),
                              rcGroupRates ? &(*rcGroupRates) : nullptr,
                              [this, curTime](SummaryState& state)
                              {
                                  this->applyParentOnlySummaryValues_(state, curTime);
                              });
        }

        // The SummaryState is now fully populated for this step, so the
        // DUMPFLUX summary snapshot can be taken.
        this->sampleFluxSummary_(isSubStep, curTime);
    }

    //! \brief Writes the initial FIP report as configured in RPTSOL.
    void writeInitialFIPReport()
    {
        const auto& gridView = simulator_.vanguard().gridView();
        const int num_interior = detail::
            countLocalInteriorCellsGridView(gridView);

        this->outputModule_->
            allocBuffers(num_interior, 0, false, false,
                         /*forceRestartFieldAllocation=*/false);

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int dofIdx = 0; dofIdx < num_interior; ++dofIdx) {
            const auto& intQuants = simulator_.model().intensiveQuantities(dofIdx, /*timeIdx=*/0);
            const auto totVolume = simulator_.model().dofTotalVolume(dofIdx);

            this->outputModule_->updateFluidInPlace(dofIdx, intQuants, totVolume);
        }

        // We always calculate the initial fip values as it may be used by various
        // keywords in the Schedule, e.g. FIP=2 in RPTSCHED but no FIP in RPTSOL
        outputModule_->calc_initial_inplace(simulator_.gridView().comm());

        // check if RPTSOL entry has FIP output
        const auto& fip = simulator_.vanguard().eclState().getEclipseConfig().fip();
        if (fip.output(FIPConfig::OutputField::FIELD) ||
            fip.output(FIPConfig::OutputField::RESV))
        {
            OPM_TIMEBLOCK(outputFipLogAndFipresvLog);

            const auto start_time = boost::posix_time::from_time_t
                (TimeService::to_time_t(simulator_.vanguard().schedule().getStartTime()));

            if (this->collectOnIORank_.isIORank()) {
                this->inplace_ = *this->outputModule_->initialInplace();

                this->outputModule_->
                    outputFipAndResvLog(this->inplace_, 0, 0.0, start_time,
                                        false, simulator_.gridView().comm());
            }
        }

        outputModule_->outputFipAndResvLogToCSV(0, false, simulator_.gridView().comm());
    }

    void writeReports(const SimulatorTimer& timer)
    {
        if (! this->collectOnIORank_.isIORank()) {
            return;
        }

        // SimulatorTimer::reportStepNum() is the simulator's zero-based
        // "episode index".  This is generally the index value needed to
        // look up objects in the Schedule container.  That said, function
        // writeReports() is invoked at the *beginning* of a report
        // step/episode which means we typically need the objects from the
        // *previous* report step/episode.  We therefore need special case
        // handling for reportStepNum() == 0 in base runs and
        // reportStepNum() <= restart step in restarted runs.
        const auto firstStep = this->initialStep();
        const auto simStep =
            std::max(timer.reportStepNum() - 1, firstStep);

        const auto& rpt = this->schedule_[simStep].rpt_config();

        if (rpt.contains("WELSPECS") && (rpt.at("WELSPECS") > 0)) {
            // Requesting a well specification report is valid at all times,
            // including reportStepNum() == initialStep().
            this->writeWellspecReport(timer);
        }

        if (timer.reportStepNum() == firstStep) {
            // No dynamic flows at the beginning of the initialStep().
            return;
        }

        if (rpt.contains("WELLS") && rpt.at("WELLS") > 0) {
            this->writeWellflowReport(timer, simStep, rpt.at("WELLS"));
        }

        this->outputModule_->outputFipAndResvLog(this->inplace_,
                                                 timer.reportStepNum(),
                                                 timer.simulationTimeElapsed(),
                                                 timer.currentDateTime(),
                                                 /* isSubstep = */ false,
                                                 simulator_.gridView().comm());

        OpmLog::note("");   // Blank line after all reports.
    }

    void writeOutput(data::Solution&& localCellData, const bool isSubStep, const bool isForcedFinalOutput)
    {
        OPM_TIMEBLOCK(writeOutput);

        const int reportStepNum = simulator_.episodeIndex() + 1;

        this->prepareLocalCellData(isSubStep, reportStepNum);
        this->outputModule_->outputErrorLog(simulator_.gridView().comm());

        // output using eclWriter if enabled
        auto localWellData = simulator_.problem().wellModel().wellData();
        auto localGroupAndNetworkData = simulator_.problem().wellModel()
            .groupAndNetworkData(reportStepNum);

        auto localAquiferData = simulator_.problem().aquiferModel().aquiferData();
        auto localWellTestState = simulator_.problem().wellModel().wellTestState();

        const bool isFlowsn = this->outputModule_->getFlows().hasFlowsn();
        auto flowsn = this->outputModule_->getFlows().getFlowsn();

        const bool isFloresn = this->outputModule_->getFlows().hasFloresn();
        auto floresn = this->outputModule_->getFlows().getFloresn();

        // assignToSolution() moves the FLORES buffers into the restart solution,
        // so DUMPFLUX has to sample its boundary faces before that happens.
        if (! isSubStep) {
            this->captureFluxDumperFlores_();
        }

        if (! isSubStep || Parameters::Get<Parameters::EnableWriteAllSolutions>()) {

            if (localCellData.empty()) {
                this->outputModule_->assignToSolution(localCellData);
            }

            // Add cell data to perforations for RFT output
            this->outputModule_->addRftDataToWells(localWellData,
                                                   reportStepNum,
                                                   simulator_.gridView().comm());
        }

        if (this->collectOnIORank_.isParallel() ||
            this->collectOnIORank_.doesNeedReordering())
        {
            // Note: We don't need WBP (well-block averaged pressures) or
            // inter-region flow rate values in order to create restart file
            // output.  There's consequently no need to collect those
            // properties on the I/O rank.

            this->collectOnIORank_.collect(localCellData,
                                           this->outputModule_->getBlockData(),
                                           this->outputModule_->getExtraBlockData(),
                                           localWellData,
                                           /* wbpData = */ {},
                                           localGroupAndNetworkData,
                                           localAquiferData,
                                           localWellTestState,
                                           /* interRegFlows = */ {},
                                           flowsn,
                                           floresn,
                                           /* lgrBlockData = */ {});
            if (this->collectOnIORank_.isIORank()) {
                this->outputModule_->assignGlobalFieldsToSolution(this->collectOnIORank_.globalCellData());
            }
        } else {
            this->outputModule_->assignGlobalFieldsToSolution(localCellData);
        }

        // Every rank takes part: a rank can only report on the cells it holds,
        // so forming a boundary record is collective even though only the IO
        // rank keeps one. Placed after the collect() above, which is what makes
        // the gathered NNC fluxes available to it.
        this->updateFluxDumpers_(reportStepNum, isSubStep);

        if (this->collectOnIORank_.isIORank()) {
            const Scalar curTime = simulator_.time() + simulator_.timeStepSize();
            const Scalar nextStepSize = simulator_.problem().nextTimeStepSize();
            std::optional<int> timeStepIdx;
            if (Parameters::Get<Parameters::EnableWriteAllSolutions>()) {
                timeStepIdx = simulator_.timeStepIndex();
            }
            this->doWriteOutput(reportStepNum, timeStepIdx, isSubStep,
                                isForcedFinalOutput,
                                std::move(localCellData),
                                std::move(localWellData),
                                std::move(localGroupAndNetworkData),
                                std::move(localAquiferData),
                                std::move(localWellTestState),
                                this->actionState(),
                                this->udqState(),
                                this->summaryState(),
                                this->simulator_.problem().thresholdPressure().getRestartVector(),
                                curTime, nextStepSize,
                                Parameters::Get<Parameters::EclOutputDoublePrecision>(),
                                isFlowsn, std::move(flowsn),
                                isFloresn, std::move(floresn));
        }
    }

    //! \brief Warn once if the boundary component masses disagree with the
    //!        residual the simulator assembled for the same faces.
    void reportFluxMassMismatch_(const std::vector<double>& written,
                                 const std::vector<double>& reference)
    {
        if (written.size() != reference.size()) {
            return;
        }

        if (this->fluxMassMismatchReported_) {
            return;
        }

        double worst = 0.0;
        double scale = 0.0;
        for (std::size_t i = 0; i < written.size(); ++i) {
            worst = std::max(worst, std::abs(written[i] - reference[i]));
            scale = std::max(scale, std::abs(reference[i]));
        }

        if (!(scale > 0.0)) {
            return;
        }

        this->fluxMassMismatchReported_ = true;

        if (worst > 1.0e-8 * scale) {
            OpmLog::warning(fmt::format("DUMPFLUX: boundary component masses differ from the "
                                        "assembled residual by up to {:.6e} against a largest "
                                        "value of {:.6e}, a relative {:.3e}. The stored rates "
                                        "do not describe what was transported.",
                                        worst, scale, worst / scale));
        }
        else {
            OpmLog::note(fmt::format("DUMPFLUX: boundary component masses agree with the "
                                     "assembled residual to a relative {:.3e}.",
                                     worst / scale));
        }
    }

    //! \brief Sample the DUMPFLUX boundary rates for the time step that just finished.
    //! \details The FLUX payload is declared as averaged over a report step. The
    //!          values are taken straight from the linearizer, which holds the
    //!          converged fluxes of the step and - unlike the output module's
    //!          buffers - is not tied to the restart output cadence.
    void sampleFluxDumperRates(const Scalar dt)
    {
        if (this->fluxDumpers_.empty()) {
            return;
        }

        if (!(dt > Scalar{0})) {
            return;
        }

        const auto& floresInfo = this->simulator_.problem().model().linearizer().getFloresInfo();
        if (floresInfo.empty()) {
            return;
        }

        const auto& flowsInfo = this->simulator_.problem().model().linearizer().getFlowsInfo();

        const auto& vanguard = this->simulator_.vanguard();

        auto floresValue = [&floresInfo, &vanguard](const int globalCell,
                                                    const FaceDir::DirEnum dir,
                                                    const int eqIdx) -> double
        {
            const auto cell = vanguard.compressedIndex(globalCell);
            if (cell < 0
                || static_cast<std::size_t>(cell) >= static_cast<std::size_t>(floresInfo.size()))
            {
                return 0.0;
            }

            const auto faceId = FaceDir::ToIntersectionIndex(dir);
            for (const auto& info : floresInfo[cell]) {
                if (info.faceId == faceId) {
                    return info.flow[eqIdx];
                }
            }

            return 0.0;
        };

        // The component surface-volume flux the simulator actually put into the
        // residual for this face, i.e. the quantity the boundary rate has to
        // reproduce. FLORES carries the phase volumetric flux, FLOWS the
        // component flux after the Rs/Rv split.
        auto flowsValue = [&flowsInfo, &vanguard](const int globalCell,
                                                  const FaceDir::DirEnum dir,
                                                  const int eqIdx) -> double
        {
            const auto cell = vanguard.compressedIndex(globalCell);
            if (cell < 0
                || static_cast<std::size_t>(cell) >= static_cast<std::size_t>(flowsInfo.size()))
            {
                return 0.0;
            }

            const auto faceId = FaceDir::ToIntersectionIndex(dir);
            for (const auto& info : flowsInfo[cell]) {
                if (info.faceId == faceId) {
                    return info.flow[eqIdx];
                }
            }

            return 0.0;
        };

        if (this->fluxMassSnapshots_.size() != this->fluxDumpers_.size()) {
            this->fluxMassSnapshots_.assign(this->fluxDumpers_.size(), {});
        }

        // Black-oil state of a cell, used to turn a phase volumetric flux into
        // component masses. Returns false when the cell is not on this rank.
        const auto orientedFaceValue = [this](const auto& accessor,
                                              const FluxRegions::BoundaryFace& face,
                                              const EclIO::FluxFile::Phase phase) -> double
        {
            if (face.isNnc || face.direction == FaceDir::Unknown) {
                return 0.0;
            }

            // The face is this rank's to report only if it owns the interior
            // cell. Every rank that merely holds a copy leaves it at zero, so
            // that the sum taken when the record is formed counts it once.
            const auto& vg = this->simulator_.vanguard();
            if (!this->fluxCellOwned_(vg.compressedIndex(face.interiorGlobalCell))) {
                return 0.0;
            }

            const auto eqIdx = this->fluxEquationIndex_(phase);
            if (eqIdx < 0) {
                return 0.0;
            }

            // Values are stored per cell for the positive face directions only
            // and are oriented along the positive axis, while the FLUX file is
            // positive into the sector.
            switch (face.direction) {
            case FaceDir::XPlus:
            case FaceDir::YPlus:
            case FaceDir::ZPlus:
                return -accessor(face.interiorGlobalCell, face.direction, eqIdx);

            case FaceDir::XMinus:
                return accessor(face.exteriorGlobalCell, FaceDir::XPlus, eqIdx);

            case FaceDir::YMinus:
                return accessor(face.exteriorGlobalCell, FaceDir::YPlus, eqIdx);

            case FaceDir::ZMinus:
                return accessor(face.exteriorGlobalCell, FaceDir::ZPlus, eqIdx);

            default:
                return 0.0;
            }
        };

        for (std::size_t i = 0; i < this->fluxDumpers_.size(); ++i) {
            // Convert to COMPONENT mass rates here, where the state of the cell
            // the flow actually comes from is known.
            auto massSnapshot = this->fluxDumpers_[i].makeFaceMajorRates(
                [&orientedFaceValue, &floresValue, this]
                (const FluxRegions::BoundaryFace& face, const EclIO::FluxFile::Phase component)
                {
                    return this->fluxComponentMass_(
                        face, component,
                        [&orientedFaceValue, &floresValue, &face](const EclIO::FluxFile::Phase p)
                        {
                            return orientedFaceValue(floresValue, face, p);
                        });
                });

            if (std::any_of(massSnapshot.begin(), massSnapshot.end(),
                            [](const double v) { return v != 0.0; }))
            {
                this->fluxMassUsable_ = true;
            }

            // The same quantity read straight out of the residual, where the
            // simulator has already done the Rs/Rv split. Any disagreement
            // means the reconstruction above does not describe what was
            // actually transported.
            if (!flowsInfo.empty()) {
                auto reference = this->fluxDumpers_[i].makeFaceMajorRates(
                    [&orientedFaceValue, &flowsValue, &vanguard, this]
                    (const FluxRegions::BoundaryFace& face, const EclIO::FluxFile::Phase component)
                    {
                        const auto interiorCell = vanguard.compressedIndex(face.interiorGlobalCell);
                        if (interiorCell < 0) {
                            return 0.0;
                        }

                        const auto pvtRegionIdx = static_cast<unsigned>(
                            this->simulator_.problem().pvtRegionIndex(interiorCell));

                        const auto phaseIdx = (component == EclIO::FluxFile::Phase::Oil)
                            ? FluidSystem::oilPhaseIdx
                            : ((component == EclIO::FluxFile::Phase::Gas)
                               ? FluidSystem::gasPhaseIdx
                               : FluidSystem::waterPhaseIdx);

                        return orientedFaceValue(flowsValue, face, component)
                            * FluidSystem::referenceDensity(phaseIdx, pvtRegionIdx);
                    });

                this->reportFluxMassMismatch_(massSnapshot, reference);
            }

            this->fluxMassSnapshots_[i].push_back(std::move(massSnapshot));
        }

        this->fluxRateTimeWeights_.push_back(static_cast<double>(dt));
    }

    void beginRestart()
    {
        const auto enablePCHysteresis = simulator_.problem().materialLawManager()->enablePCHysteresis();
        const auto enableNonWettingHysteresis = simulator_.problem().materialLawManager()->enableNonWettingHysteresis();
        const auto enableWettingHysteresis = simulator_.problem().materialLawManager()->enableWettingHysteresis();
        const auto oilActive = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
        const auto gasActive = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
        const auto waterActive = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
        const auto enableSwatinit = simulator_.vanguard().eclState().fieldProps().has_double("SWATINIT");

        std::vector<RestartKey> solutionKeys {
            {"PRESSURE", UnitSystem::measure::pressure},
            {"SWAT",     UnitSystem::measure::identity,    waterActive},
            {"SGAS",     UnitSystem::measure::identity,    gasActive},
            {"TEMP",     UnitSystem::measure::temperature, enableEnergy},
            {"SSOLVENT", UnitSystem::measure::identity,    enableSolvent},

            {"RS",  UnitSystem::measure::gas_oil_ratio, FluidSystem::enableDissolvedGas()},
            {"RV",  UnitSystem::measure::oil_gas_ratio, FluidSystem::enableVaporizedOil()},
            {"RVW", UnitSystem::measure::oil_gas_ratio, FluidSystem::enableVaporizedWater()},
            {"RSW", UnitSystem::measure::gas_oil_ratio, FluidSystem::enableDissolvedGasInWater()},

            {"SGMAX", UnitSystem::measure::identity, enableNonWettingHysteresis && oilActive && gasActive},
            {"SHMAX", UnitSystem::measure::identity, enableWettingHysteresis && oilActive && gasActive},

            {"SOMAX", UnitSystem::measure::identity,
             (enableNonWettingHysteresis && oilActive && waterActive)
             || simulator_.problem().vapparsActive(simulator_.episodeIndex())},

            {"SOMIN", UnitSystem::measure::identity, enablePCHysteresis && oilActive && gasActive},
            {"SWHY1", UnitSystem::measure::identity, enablePCHysteresis && oilActive && waterActive},
            {"SWMAX", UnitSystem::measure::identity, enableWettingHysteresis && oilActive && waterActive},

            {"PPCW", UnitSystem::measure::pressure, enableSwatinit},
        };

        {
            const auto& tracers = simulator_.vanguard().eclState().tracer();

            for (const auto& tracer : tracers) {
                const auto enableSolTracer =
                    ((tracer.phase == Phase::GAS) && FluidSystem::enableDissolvedGas()) ||
                    ((tracer.phase == Phase::OIL) && FluidSystem::enableVaporizedOil());

                solutionKeys.emplace_back(tracer.fname(), UnitSystem::measure::identity, true);
                solutionKeys.emplace_back(tracer.sname(), UnitSystem::measure::identity, enableSolTracer);
            }
        }

        const auto& inputThpres = eclState().getSimulationConfig().getThresholdPressure();
        const std::vector<RestartKey> extraKeys {
            {"OPMEXTRA", UnitSystem::measure::identity, false},
            {"THRESHPR", UnitSystem::measure::pressure, inputThpres.active()},
        };

        const auto& gridView = this->simulator_.vanguard().gridView();
        const auto numElements = gridView.size(/*codim=*/0);

        // Try to load restart step 0 to calculate initial FIP
        {
            this->outputModule_->allocBuffers(numElements,
                                              0,
                                              /*isSubStep = */false,
                                              /*log = */      false,
                                              /*forceRestartFieldAllocation = */true);

            const auto restartSolution =
                loadParallelRestartSolution(this->eclIO_.get(),
                                            solutionKeys, gridView.comm(), 0);

            if (!restartSolution.empty()) {
                for (auto elemIdx = 0*numElements; elemIdx < numElements; ++elemIdx) {
                    const auto globalIdx = this->collectOnIORank_.localIdxToGlobalIdx(elemIdx);
                    this->outputModule_->setRestart(restartSolution, elemIdx, globalIdx);
                }

                this->simulator_.problem().readSolutionFromOutputModule(0, true);
                this->simulator_.problem().temperatureModel().init();
                ElementContext elemCtx(this->simulator_);
                for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                    elemCtx.updatePrimaryStencil(elem);
                    elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                    this->outputModule_->updateFluidInPlace(elemCtx);
                }

                this->outputModule_->calc_initial_inplace(this->simulator_.gridView().comm());
            }
        }

        {
            // The episodeIndex is rewound one step back before calling
            // beginRestart() and cannot be used here.  We just ask the
            // initconfig directly to be sure that we use the correct index.
            const auto restartStepIdx = this->simulator_.vanguard()
                .eclState().getInitConfig().getRestartStep();

            this->outputModule_->allocBuffers(numElements,
                                              restartStepIdx,
                                              /*isSubStep = */false,
                                              /*log = */      false,
                                              /*forceRestartFieldAllocation = */true);
        }

        {
            const auto restartValues =
                loadParallelRestart(this->eclIO_.get(),
                                    this->actionState(),
                                    this->summaryState(),
                                    solutionKeys, extraKeys, gridView.comm());

            for (auto elemIdx = 0*numElements; elemIdx < numElements; ++elemIdx) {
                const auto globalIdx = this->collectOnIORank_.localIdxToGlobalIdx(elemIdx);
                this->outputModule_->setRestart(restartValues.solution, elemIdx, globalIdx);
            }

            auto& tracer_model = simulator_.problem().tracerModel();
            for (int tracer_index = 0; tracer_index < tracer_model.numTracers(); ++tracer_index) {
                // Free tracers
                {
                    const auto& free_tracer_name = tracer_model.fname(tracer_index);
                    const auto& free_tracer_solution = restartValues.solution
                        .template data<double>(free_tracer_name);

                    for (auto elemIdx = 0*numElements; elemIdx < numElements; ++elemIdx) {
                        const auto globalIdx = this->collectOnIORank_.localIdxToGlobalIdx(elemIdx);
                        tracer_model.setFreeTracerConcentration
                            (tracer_index, elemIdx, free_tracer_solution[globalIdx]);
                    }
                }

                // Solution tracer (only if DISGAS/VAPOIL are active for gas/oil tracers)
                if ((tracer_model.phase(tracer_index) == Phase::GAS && FluidSystem::enableDissolvedGas()) ||
                    (tracer_model.phase(tracer_index) == Phase::OIL && FluidSystem::enableVaporizedOil()))
                {
                    tracer_model.setEnableSolTracers(tracer_index, true);

                    const auto& sol_tracer_name = tracer_model.sname(tracer_index);
                    const auto& sol_tracer_solution = restartValues.solution
                        .template data<double>(sol_tracer_name);

                    for (auto elemIdx = 0*numElements; elemIdx < numElements; ++elemIdx) {
                        const auto globalIdx = this->collectOnIORank_.localIdxToGlobalIdx(elemIdx);
                        tracer_model.setSolTracerConcentration
                            (tracer_index, elemIdx, sol_tracer_solution[globalIdx]);
                    }
                }
                else {
                    tracer_model.setEnableSolTracers(tracer_index, false);

                    for (auto elemIdx = 0*numElements; elemIdx < numElements; ++elemIdx) {
                        tracer_model.setSolTracerConcentration(tracer_index, elemIdx, 0.0);
                    }
                }
            }

            if (inputThpres.active()) {
                const_cast<Simulator&>(this->simulator_)
                    .problem().thresholdPressure()
                    .setFromRestart(restartValues.getExtra("THRESHPR"));
            }

            restartTimeStepSize_ = restartValues.getExtra("OPMEXTRA")[0];
            if (restartTimeStepSize_ <= 0) {
                restartTimeStepSize_ = std::numeric_limits<double>::max();
            }

            // Initialize the well model from restart values
            this->simulator_.problem().wellModel()
                .initFromRestartFile(restartValues);

            if (!restartValues.aquifer.empty()) {
                this->simulator_.problem().mutableAquiferModel()
                    .initFromRestart(restartValues.aquifer);
            }
        }
    }

    void endRestart()
    {
        // Calculate initial in-place volumes.
        // Does nothing if they have already been calculated,
        // e.g. from restart data at T=0.
        this->outputModule_->calc_initial_inplace(this->simulator_.gridView().comm());

        if (this->collectOnIORank_.isIORank()) {
            if (const auto* iip = this->outputModule_->initialInplace(); iip != nullptr) {
                this->inplace_ = *iip;
            }
        }
    }

    const OutputModule& outputModule() const
    { return *outputModule_; }

    OutputModule& mutableOutputModule() const
    { return *outputModule_; }

    Scalar restartTimeStepSize() const
    { return restartTimeStepSize_; }

    template <class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(*outputModule_);
    }

private:
    static bool enableEclOutput_()
    {
        static bool enable = Parameters::Get<Parameters::EnableEclOutput>();
        return enable;
    }

    const EclipseState& eclState() const
    { return simulator_.vanguard().eclState(); }

    SummaryState& summaryState()
    { return simulator_.vanguard().summaryState(); }

    const SummaryState& summaryState() const
    { return simulator_.vanguard().summaryState(); }

    Action::State& actionState()
    { return simulator_.vanguard().actionState(); }

    UDQState& udqState()
    { return simulator_.vanguard().udqState(); }

    const Schedule& schedule() const
    { return simulator_.vanguard().schedule(); }

    /// Collect the reservoir coupling data Summary::eval() needs: on a
    /// master, the slaves' group rates; on a slave, the injection targets and
    /// production limits in force for its slave groups.  Returns nullopt for
    /// non-RC simulations.
    std::optional<data::ReservoirCouplingGroupRates> collectReservoirCouplingGroupRates_()
    {
#ifdef RESERVOIR_COUPLING_ENABLED
        // Guard: only BlackoilWellModel has reservoir coupling support.
        // CompWellModel (compositional) does not, so we use if constexpr
        // to avoid compilation errors when EclWriter is instantiated with
        // a compositional TypeTag.
        using WellModelType = std::remove_cvref_t<
            decltype(simulator_.problem().wellModel())>;
        if constexpr (requires(WellModelType& wm) { wm.isReservoirCouplingMaster(); }) {
            auto& wellModel = simulator_.problem().wellModel();
            if (wellModel.isReservoirCouplingMaster()) {
                return wellModel.reservoirCouplingMaster()
                    .collectGroupRatesForSummary();
            }
            if (wellModel.isReservoirCouplingSlave()) {
                auto rates = data::ReservoirCouplingGroupRates{};
                for (const auto& [group, targets] :
                         wellModel.reservoirCouplingSlave().effectiveInjectionTargets())
                {
                    for (const auto& [phase, target] : targets) {
                        rates.injection_targets[group][phase] = static_cast<double>(target);
                    }
                }
                for (const auto& [group, limits] :
                         wellModel.reservoirCouplingSlave().effectiveProductionTargets())
                {
                    for (const auto& [cmode, limit] : limits) {
                        rates.production_targets[group][cmode] = static_cast<double>(limit);
                    }
                }
                return rates;
            }
        }
#endif
        return std::nullopt;
    }

    void prepareLocalCellData(const bool isSubStep,
                              const int  reportStepNum)
    {
        OPM_TIMEBLOCK(prepareLocalCellData);

        if (this->outputModule_->localDataValid()) {
            return;
        }

        const auto& gridView = simulator_.vanguard().gridView();
        const bool log = this->collectOnIORank_.isIORank();

        const int num_interior = detail::
            countLocalInteriorCellsGridView(gridView);
        const bool writeAllSolutions =
            Parameters::Get<Parameters::EnableWriteAllSolutions>();

        // EclipseIO writes restart output for every positive time-step index in
        // write-all mode, independently of the schedule's BASIC/FREQ settings.
        const bool forceRestartFieldAllocation =
            writeAllSolutions && (simulator_.timeStepIndex() > 0);
        this->outputModule_->
            allocBuffers(num_interior, reportStepNum,
                         isSubStep && !writeAllSolutions,
                         log, forceRestartFieldAllocation);

        ElementContext elemCtx(simulator_);

        OPM_BEGIN_PARALLEL_TRY_CATCH();

        {
            OPM_TIMEBLOCK(prepareCellBasedData);

            this->outputModule_->prepareDensityAccumulation();
            this->outputModule_->setupExtractors(isSubStep, reportStepNum);
            for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                this->outputModule_->processElement(elemCtx);
                this->outputModule_->processElementBlockData(elemCtx);
            }
            this->outputModule_->clearExtractors();
        }

        {
            OPM_TIMEBLOCK(prepareFluidInPlace);

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int dofIdx = 0; dofIdx < num_interior; ++dofIdx) {
                const auto& intQuants = simulator_.model().intensiveQuantities(dofIdx, /*timeIdx=*/0);
                const auto totVolume = simulator_.model().dofTotalVolume(dofIdx);

                this->outputModule_->updateFluidInPlace(dofIdx, intQuants, totVolume);
            }
        }

        OPM_END_PARALLEL_TRY_CATCH("EclWriter::prepareLocalCellData() failed: ",
                                   this->simulator_.vanguard().grid().comm());

        // Propagate rank-local exceptions before entering output collectives.
        this->outputModule_->accumulateDensityParallel();
        this->outputModule_->validateLocalData();
    }

    void captureLocalFluxData()
    {
        OPM_TIMEBLOCK(captureLocalData);

        const auto& gridView = this->simulator_.vanguard().gridView();
        const auto timeIdx = 0u;

        auto elemCtx = ElementContext { this->simulator_ };

        const auto elemMapper = ElementMapper { gridView, Dune::mcmgElementLayout() };
        const auto activeIndex = [&elemMapper](const Element& e)
        {
            return elemMapper.index(e);
        };

        const auto cartesianIndex = [this](const int elemIndex)
        {
            return this->cartMapper_.cartesianIndex(elemIndex);
        };

        this->outputModule_->initializeFluxData();

        OPM_BEGIN_PARALLEL_TRY_CATCH();

        for (const auto& elem : elements(gridView, Dune::Partitions::interiorBorder)) {
            elemCtx.updateStencil(elem);
            elemCtx.updateIntensiveQuantities(timeIdx);
            elemCtx.updateExtensiveQuantities(timeIdx);

            this->outputModule_->processFluxes(elemCtx, activeIndex, cartesianIndex);
        }

        OPM_END_PARALLEL_TRY_CATCH("EclWriter::captureLocalFluxData() failed: ",
                                   this->simulator_.vanguard().grid().comm())

        this->outputModule_->finalizeFluxData();
    }

    void writeWellspecReport(const SimulatorTimer& timer) const
    {
        const auto changedWells = this->schedule_
            .changed_wells(timer.reportStepNum(), this->initialStep());

        const auto changedWellLists = this->schedule_
            .changedWellLists(timer.reportStepNum(), this->initialStep());

        if (changedWells.empty() && !changedWellLists) {
            return;
        }

        this->outputModule_->outputWellspecReport(changedWells,
                                                  changedWellLists,
                                                  timer.reportStepNum(),
                                                  timer.simulationTimeElapsed(),
                                                  timer.currentDateTime());
    }

    void writeWellflowReport(const SimulatorTimer& timer,
                             const int             simStep,
                             const int             wellsRequest) const
    {
        this->outputModule_->outputTimeStamp("WELLS",
                                             timer.simulationTimeElapsed(),
                                             timer.reportStepNum(),
                                             timer.currentDateTime());

        const auto wantConnData = wellsRequest > 1;

        this->outputModule_->outputProdLog(simStep, wantConnData);
        this->outputModule_->outputInjLog(simStep, wantConnData);
        this->outputModule_->outputCumLog(simStep, wantConnData);
        this->outputModule_->outputMSWLog(simStep);
    }

    int initialStep() const
    {
        const auto& initConfig = this->eclState().cfg().init();

        return initConfig.restartRequested()
            ? initConfig.getRestartStep()
            : 0;
    }

    int fluxPhaseMask_() const
    {
        auto mask = 0;

        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            mask |= static_cast<int>(EclIO::FluxFile::Phase::Oil);
        }
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            mask |= static_cast<int>(EclIO::FluxFile::Phase::Water);
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            mask |= static_cast<int>(EclIO::FluxFile::Phase::Gas);
        }

        return mask;
    }

    // Black-oil state of a cell, used to turn a phase volumetric flux into
    // component masses. Not filled when the cell is not on this rank.
    struct FluxCellPvt
    {
        std::array<double, 3> invB{{0.0, 0.0, 0.0}};  // canonical phase index
        double rs = 0.0;
        double rv = 0.0;
        bool ok = false;
    };

    //! \brief Whether the fluid system is black-oil, which DUMPFLUX needs:
    //!        a FLUX file holds formation volume factors, Rs and Rv, which a
    //!        compositional run does not have.
    static constexpr bool fluxBlackOil_ =
        requires { FluidSystem::referenceDensity(0u, 0u); };

    FluxCellPvt fluxCellPvt_(const int globalCell) const
    {
        FluxCellPvt out{};
        if constexpr (!fluxBlackOil_) {
            return out;
        }
        else {
            const auto cell = this->simulator_.vanguard().compressedIndex(globalCell);
            if (cell < 0) {
                return out;
            }

            const auto* intQuants =
                this->simulator_.model().cachedIntensiveQuantities(cell, /*timeIdx=*/0);
            if (intQuants == nullptr) {
                return out;
            }

            const auto& fs = intQuants->fluidState();
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (FluidSystem::phaseIsActive(phaseIdx)) {
                    out.invB[phaseIdx] = getValue(fs.invB(phaseIdx));
                }
            }
            if (FluidSystem::enableDissolvedGas()) {
                out.rs = getValue(fs.Rs());
            }
            if (FluidSystem::enableVaporizedOil()) {
                out.rv = getValue(fs.Rv());
            }

            out.ok = true;
            return out;
        }
    }

    // Convert the phase volumetric fluxes across one boundary face, supplied by
    // rateOf() and positive into the sector, into the mass of one component.
    //
    // This has to happen where the state of the cell the flow comes from is
    // known. A reduced run cannot do it for inflow, because that upstream cell
    // lies outside its grid.
    //
    // The split matters for live oil and wet gas: the oil phase carries gas
    // component mass through Rs and the gas phase carries oil component mass
    // through Rv, so multiplying a phase flux by the phase density and calling
    // the result "oil" would put a sizeable part of it into the wrong
    // conservation equation.
    //
    // The reference densities are taken from the INTERIOR cell, which is the
    // one the consumer will use when it converts these masses back to surface
    // volumes, so the round trip is exact even where the two sides of the face
    // are in different PVT regions.
    template <typename RateOf>
    double fluxComponentMass_(const FluxRegions::BoundaryFace& face,
                              const EclIO::FluxFile::Phase component,
                              RateOf&& rateOf) const
    {
        const auto interiorCell =
            this->simulator_.vanguard().compressedIndex(face.interiorGlobalCell);
        if (interiorCell < 0) {
            return 0.0;
        }

        const auto pvtRegionIdx =
            static_cast<unsigned>(this->simulator_.problem().pvtRegionIndex(interiorCell));

        // Upwind state for one phase: positive is into the sector, so inflow
        // comes from the exterior cell.
        const auto upwindOf = [&face, this](const double volRate)
        {
            return this->fluxCellPvt_((volRate > 0.0) ? face.exteriorGlobalCell
                                                      : face.interiorGlobalCell);
        };

        const bool hasOil = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
        const bool hasGas = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
        const bool hasWat = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);

        const double qo = hasOil ? rateOf(EclIO::FluxFile::Phase::Oil) : 0.0;
        const double qg = hasGas ? rateOf(EclIO::FluxFile::Phase::Gas) : 0.0;
        const double qw = hasWat ? rateOf(EclIO::FluxFile::Phase::Water) : 0.0;

        // Surface volume flux of each phase, i.e. reservoir volume times the
        // inverse formation volume factor upstream.
        double sO = 0.0;
        double sG = 0.0;
        double sW = 0.0;
        double rs = 0.0;
        double rv = 0.0;

        if (hasOil && (qo != 0.0)) {
            const auto up = upwindOf(qo);
            if (up.ok) {
                sO = qo * up.invB[FluidSystem::oilPhaseIdx];
                rs = up.rs;
            }
        }
        if (hasGas && (qg != 0.0)) {
            const auto up = upwindOf(qg);
            if (up.ok) {
                sG = qg * up.invB[FluidSystem::gasPhaseIdx];
                rv = up.rv;
            }
        }
        if (hasWat && (qw != 0.0)) {
            const auto up = upwindOf(qw);
            if (up.ok) {
                sW = qw * up.invB[FluidSystem::waterPhaseIdx];
            }
        }

        if constexpr (!fluxBlackOil_) {
            return 0.0;
        }
        else {
            switch (component) {
            case EclIO::FluxFile::Phase::Oil:
                return (sO + rv * sG)
                    * FluidSystem::referenceDensity(FluidSystem::oilPhaseIdx, pvtRegionIdx);

            case EclIO::FluxFile::Phase::Gas:
                return (sG + rs * sO)
                    * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIdx);

            case EclIO::FluxFile::Phase::Water:
                return sW
                    * FluidSystem::referenceDensity(FluidSystem::waterPhaseIdx, pvtRegionIdx);
            }

            return 0.0;
        }
    }

    int fluxComponentIndex_(const EclIO::FluxFile::Phase phase) const
    {
        switch (phase) {
        case EclIO::FluxFile::Phase::Oil:
            return FluidSystem::oilCompIdx;
        case EclIO::FluxFile::Phase::Water:
            return FluidSystem::waterCompIdx;
        case EclIO::FluxFile::Phase::Gas:
            return FluidSystem::gasCompIdx;
        }

        return FluidSystem::oilCompIdx;
    }

    int fluxEquationIndex_(const EclIO::FluxFile::Phase phase) const
    {
        switch (phase) {
        case EclIO::FluxFile::Phase::Oil:
            return Indices::oilEnabled
                ? Indices::conti0EqIdx + FluidSystem::canonicalToActiveCompIdx(FluidSystem::oilCompIdx)
                : -1;
        case EclIO::FluxFile::Phase::Water:
            return Indices::waterEnabled
                ? Indices::conti0EqIdx + FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx)
                : -1;
        case EclIO::FluxFile::Phase::Gas:
            return Indices::gasEnabled
                ? Indices::conti0EqIdx + FluidSystem::canonicalToActiveCompIdx(FluidSystem::gasCompIdx)
                : -1;
        }

        return -1;
    }

    EclIO::FluxFile::Mode fluxOutputMode_(const std::string& fluxType) const
    {
        if (fluxType == "PRESSURE") {
            return EclIO::FluxFile::Mode::Pressure;
        }
        if (fluxType == "BOTH") {
            return EclIO::FluxFile::Mode::Both;
        }
        return EclIO::FluxFile::Mode::Flux;
    }

    static std::pair<int, int> normalizedNncPair_(const int c1, const int c2)
    {
        return (c1 <= c2)
            ? std::make_pair(c1, c2)
            : std::make_pair(c2, c1);
    }

    //! \brief Append the non-neighbour connections the grid itself makes.
    //!
    //! \details The deck's NNC keyword and PINCH are not the only sources of
    //!   connections between cells that are not axis neighbours. A fault with
    //!   throw joins each cell on one side to whichever cells face it on the
    //!   other, several layers up or down, and nothing in the input lists
    //!   those: the grid works them out from the corner-point geometry. Where
    //!   such a fault runs along a sector boundary, every one of them is a
    //!   boundary face of the sector, and leaving them out seals the fault in
    //!   the reduced run however open it is in the full one. On one field that
    //!   left a column of cells against a fault 8 bar off and drove the water
    //!   cut of a producer four cells away steadily higher.
    //!
    //!   Pairs of vertical neighbours with nothing but inactive cells between
    //!   them are connected directly and are not NNCs, so they are left out,
    //!   as are axis neighbours; both are what the Cartesian faces cover.
    //!   Duplicates of the input connections are harmless: the region builder
    //!   discards a face it has already seen.
    //!
    //!   Walks the serial I/O grid, and so only on the I/O rank.
    void appendGridNncPairs_(std::vector<int>& nncPairs) const
    {
        if ((this->equilGrid_ == nullptr) || (this->equilCartMapper_ == nullptr)) {
            return;
        }

        using GlobalGridView = typename EquilGrid::LeafGridView;
        using GlobalElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GlobalGridView>;

        const GlobalGridView gridView = this->equilGrid_->leafGridView();
        const GlobalElementMapper elemMapper { gridView, Dune::mcmgElementLayout() };
        const auto& cartMapper = *this->equilCartMapper_;
        const auto dims = cartMapper.cartesianDimensions();

        std::unordered_set<int> active;
        for (const auto& elem : elements(gridView)) {
            active.insert(cartMapper.cartesianIndex(elemMapper.index(elem)));
        }

        const auto ijkOf = [&dims](int cell)
        {
            const int i = cell % dims[0];
            cell /= dims[0];
            return std::array<int, 3>{ i, cell % dims[1], cell / dims[1] };
        };

        const auto directlyConnected = [&](const int c1, const int c2)
        {
            const auto a = ijkOf(c1);
            const auto b = ijkOf(c2);
            const int di = std::abs(a[0] - b[0]);
            const int dj = std::abs(a[1] - b[1]);
            const int dk = std::abs(a[2] - b[2]);

            if (di + dj + dk == 1) {
                return true;
            }

            if ((di == 0) && (dj == 0)) {
                const auto step = dims[0] * dims[1];
                for (int c = std::min(c1, c2) + step; c < std::max(c1, c2); c += step) {
                    if (active.count(c) != 0) {
                        return false;
                    }
                }
                return true;
            }

            return false;
        };

        std::set<std::pair<int, int>> seen;
        for (const auto& elem : elements(gridView)) {
            const int c1 = cartMapper.cartesianIndex(elemMapper.index(elem));
            for (const auto& is : intersections(gridView, elem)) {
                if (!is.neighbor()) {
                    continue;
                }

                const int c2 = cartMapper.cartesianIndex(elemMapper.index(is.outside()));
                if ((c2 <= c1) || directlyConnected(c1, c2)) {
                    continue;
                }

                if (seen.emplace(c1, c2).second) {
                    nncPairs.push_back(c1);
                    nncPairs.push_back(c2);
                }
            }
        }
    }

    //! \brief Where each NNC sits in the writer's NNC list.
    //!
    //! \details FLORESN is indexed by that list, which holds the deck's
    //!   connections followed by the ones the grid generated, so this is the
    //!   only index that finds an NNC's flux. It is complete only once the
    //!   transmissibilities have been exported, well after the regions are set
    //!   up, hence the deferred construction.
    const std::map<std::pair<int, int>, int>& fluxNncPairToIndex_() const
    {
        if (this->fluxNncIndex_.empty()) {
            const auto& outputNnc = this->getOutputNnc();
            if (!outputNnc.empty()) {
                const auto& level0 = outputNnc.front();
                for (std::size_t n = 0; n < level0.size(); ++n) {
                    this->fluxNncIndex_.try_emplace(
                        normalizedNncPair_(static_cast<int>(level0[n].cell1),
                                           static_cast<int>(level0[n].cell2)),
                        static_cast<int>(n));
                }
            }
        }

        return this->fluxNncIndex_;
    }

    // A sector boundary cuts the grid, and a well with completions on both
    // sides of it cannot be reproduced by a reduced run: that run sees only the
    // connections inside its own region but solves the well as though it were
    // whole, so both its rates and its bottom hole pressure come out wrong.
    // Refuse to write a FLUX file that could only be used incorrectly.
    //
    // The ground outside every region counts as a region of its own here. A
    // well reaching out of the sector is just as badly served by a reduced run
    // as one reaching into a neighbouring sector.
    //
    // Only the I/O rank holds the global field properties, so only it can tell,
    // but the verdict is shared so that every rank leaves through the same door
    // and MPI shuts down cleanly.
    void checkWellsWithinSingleFluxRegion_() const
    {
        auto excType = ExceptionType::NONE;
        auto message = std::string
            {"DUMPFLUX requires every well to be completed within a single FLUXNUM "
             "region. The offending wells are listed on the I/O rank."};

        if (this->collectOnIORank_.isIORank()) {
            const auto& state = this->eclState();
            const auto& fieldProps = state.globalFieldProps();

            if (state.getIOConfig().getDumpFlux() && fieldProps.has_int("FLUXNUM")) {
                const auto offenders =
                    this->wellsStraddlingFluxRegions_(fieldProps.get_global_int("FLUXNUM"),
                                                      fieldProps.actnumRaw(),
                                                      state.gridDims().getNXYZ());

                if (!offenders.empty()) {
                    excType = ExceptionType::INVALID_ARGUMENT;
                    message = fmt::format(
                        "DUMPFLUX requires every well to be completed within a single "
                        "FLUXNUM region. A reduced run covering one region would see only "
                        "the part of a straddling well that falls inside it, but would "
                        "solve that well as though it were whole. Offending wells:\n{}\n"
                        "Either move the region boundary clear of these completions or "
                        "keep the wells out of the sector.",
                        fmt::join(offenders, "\n"));

                    // Put it in the print file as well, so the list survives
                    // the run rather than only reaching the terminal.
                    OpmLog::error(message);
                }
            }
        }

        checkForExceptionsAndThrow(excType, message,
                                   this->simulator_.vanguard().grid().comm());
    }

    std::vector<std::string>
    wellsStraddlingFluxRegions_(const std::vector<int>& regionValues,
                                const std::vector<int>& actnum,
                                const std::array<int, 3>& dims) const
    {
        const auto numCells = static_cast<std::size_t>(dims[0])
            * static_cast<std::size_t>(dims[1])
            * static_cast<std::size_t>(dims[2]);

        // Completions are added as the schedule advances, so the whole of it
        // has to be walked rather than just the final state.
        std::map<std::string, std::set<int>> wellRegions;

        const auto& sched = this->schedule();
        for (std::size_t step = 0; step < sched.size(); ++step) {
            for (const auto& well : sched.getWells(step)) {
                auto& regions = wellRegions[well.name()];

                for (const auto& conn : well.getConnections()) {
                    const auto cell = conn.global_index();
                    if (cell >= numCells) {
                        continue;
                    }

                    // A connection in a cell the grid does not have is not a
                    // connection at all, and counting the region of such a cell
                    // would report wells that are in fact perfectly placed.
                    if (!actnum.empty() && (actnum[cell] == 0)) {
                        continue;
                    }

                    regions.insert(regionValues[cell]);
                }
            }
        }

        std::vector<std::string> offenders;
        for (const auto& [name, regions] : wellRegions) {
            if (regions.size() < 2) {
                continue;
            }

            std::vector<std::string> labels;
            labels.reserve(regions.size());
            for (const auto region : regions) {
                labels.push_back((region == 0)
                                 ? std::string{"0 (outside every region)"}
                                 : std::to_string(region));
            }

            offenders.push_back(fmt::format("  {} is completed in flux regions {}",
                                            name, fmt::join(labels, ", ")));
        }

        return offenders;
    }

    //! \brief Set up one dumper per FLUXNUM region, on every rank.
    //!
    //! \details A rank can only report on the cells it holds, so every rank
    //!   has to take part in filling a boundary record even though only the IO
    //!   rank keeps one. Leaving the dumpers on the IO rank alone left it
    //!   writing zero for every face outside its own partition.
    //!
    //!   The region is a pure function of the deck, so each rank could work it
    //!   out for itself were the inputs at hand. They are not: FLUXNUM and
    //!   ACTNUM live in the global field properties, which only the root rank
    //!   may read. They are broadcast instead, and every rank then runs the
    //!   same extraction and arrives at the same region.
    void initializeFluxDumpers_()
    {
        const auto& state = this->eclState();
        const auto& io = state.getIOConfig();
        const auto& comm = this->simulator_.vanguard().grid().comm();
        const bool isRoot = this->collectOnIORank_.isIORank();

        // FLUXNUM is an ordinary region array, which a USEFLUX consumer and
        // MULTREGT among others use as well, so only DUMPFLUX makes this a
        // producer run.
        if (!io.getDumpFlux()) {
            return;
        }

        const auto dims = state.gridDims().getNXYZ();
        const auto numGlobalCells =
            static_cast<int>(static_cast<std::size_t>(dims[0]) * dims[1] * dims[2]);

        // FLUXNUM is dimensioned over every cell of the grid, active or not,
        // so it has to be read with the global accessor rather than the one
        // that returns a value per active cell. ACTNUM goes alongside it so
        // that inactive cells are kept out of the region.
        std::vector<int> regionValues;
        std::vector<int> actnum;
        std::vector<int> nncPairs;
        int haveFluxnum = 0;

        if (isRoot) {
            const auto& fieldProps = state.globalFieldProps();
            haveFluxnum = fieldProps.has_int("FLUXNUM") ? 1 : 0;

            if (haveFluxnum != 0) {
                regionValues = fieldProps.get_global_int("FLUXNUM");
                actnum = fieldProps.actnumRaw();

                const auto appendNnc = [&nncPairs](const auto& connections)
                {
                    for (const auto& nnc : connections) {
                        nncPairs.push_back(static_cast<int>(nnc.cell1));
                        nncPairs.push_back(static_cast<int>(nnc.cell2));
                    }
                };

                if (state.hasInputNNC()) {
                    appendNnc(state.getInputNNC().input());
                }
                if (state.hasPinchNNC()) {
                    appendNnc(state.getPinchNNC());
                }

                this->appendGridNncPairs_(nncPairs);
            }
        }

        comm.broadcast(&haveFluxnum, 1, 0);
        if (haveFluxnum == 0) {
            return;
        }

        const auto share = [&comm, isRoot](std::vector<int>& values)
        {
            int size = isRoot ? static_cast<int>(values.size()) : 0;
            comm.broadcast(&size, 1, 0);
            if (!isRoot) {
                values.assign(static_cast<std::size_t>(size), 0);
            }
            if (size > 0) {
                comm.broadcast(values.data(), size, 0);
            }
        };

        share(regionValues);
        share(actnum);
        share(nncPairs);

        if (static_cast<int>(regionValues.size()) != numGlobalCells) {
            return;
        }

        std::vector<std::array<int, 2>> nncConnections;
        nncConnections.reserve(nncPairs.size() / 2);
        for (std::size_t n = 0; n + 1 < nncPairs.size(); n += 2) {
            nncConnections.push_back({nncPairs[n], nncPairs[n + 1]});
        }

        const auto regions = FluxRegions::extract(dims, regionValues, actnum, nncConnections);
        if (regions.empty()) {
            return;
        }

        // The boundary state in a FLUX file is black-oil (formation volume
        // factors, Rs, Rv). Refuse rather than write a file a consumer would
        // read as zero fluxes.
        if constexpr (!fluxBlackOil_) {
            throw std::runtime_error {
                "DUMPFLUX (FLUXNUM in a producer run) is only supported "
                "for black-oil models"
            };
        }

        const auto phaseMask = this->fluxPhaseMask_();
        const auto fluxMode = this->fluxOutputMode_(io.getFluxType());

        this->fluxOutputPaths_.clear();
        this->fluxDumpers_.reserve(regions.size());
        this->fluxOutputPaths_.reserve(regions.size());

        const auto basePath = (std::filesystem::path{io.getOutputDir()} / io.getBaseName()).string();
        const auto multipleRegions = regions.size() > 1;

        int sequence = 1;
        const auto parentWells = describeFluxParentWells(this->simulator_.vanguard().schedule());
        for (const auto& region : regions) {
            this->fluxDumpers_.emplace_back(io.getBaseName(),
                                            region.regionId,
                                            dims,
                                            region,
                                            fluxMode,
                                            EclIO::FluxFile::Sampling::Averaged,
                                            phaseMask);
            this->fluxDumpers_.back().setParentWells(parentWells);

            if (multipleRegions) {
                std::ostringstream os;
                os << basePath << ".FLUX" << std::setw(4) << std::setfill('0') << sequence;
                this->fluxOutputPaths_.push_back(os.str());
            }
            else {
                this->fluxOutputPaths_.push_back(basePath + ".FLUX");
            }
            ++sequence;
        }

        OpmLog::note("DUMPFLUX bootstrap: initialized "
                     + std::to_string(this->fluxDumpers_.size())
                 + " region dumper(s) from FLUXNUM");

        this->buildFluxOwnedCells_();
        this->initializeFluxSummarySampling_();
    }

    //! \brief Mark the cells this rank is the owner of.
    //!
    //! \details A cell shows up on more than one rank, once as an interior
    //!   cell and again in the overlap of whoever borders it. A per-cell
    //!   quantity summed across ranks therefore has to be contributed by
    //!   exactly one of them, and the interior partition picks that one. Taking
    //!   a maximum instead would avoid the double count for a positive
    //!   quantity, but not for capillary pressure, which is signed.
    void buildFluxOwnedCells_()
    {
        this->fluxOwnedCell_.clear();

        if (this->fluxDumpers_.empty()) {
            return;
        }

        const auto& gridView = this->simulator_.vanguard().gridView();
        const auto& mapper = this->simulator_.model().elementMapper();

        this->fluxOwnedCell_.assign(gridView.size(/*codim=*/0), 0);
        for (const auto& elem : elements(gridView)) {
            if (elem.partitionType() == Dune::InteriorEntity) {
                this->fluxOwnedCell_[mapper.index(elem)] = 1;
            }
        }
    }

    //! \brief Whether this rank is the one that should report on a cell.
    bool fluxCellOwned_(const int compressedCell) const
    {
        return (compressedCell >= 0)
            && (static_cast<std::size_t>(compressedCell) < this->fluxOwnedCell_.size())
            && (this->fluxOwnedCell_[compressedCell] != 0);
    }

    //! \brief Gather a per-face or per-region quantity from all ranks.
    //!
    //! \details Every rank has filled the entries it owns and left the rest at
    //!   zero, so a sum collects them. Collective, hence called from the same
    //!   place on every rank.
    template <typename T>
    void fluxReduceSum_(std::vector<T>& values) const
    {
        const auto& comm = this->simulator_.vanguard().grid().comm();
        if ((comm.size() > 1) && !values.empty()) {
            comm.sum(values.data(), values.size());
        }
    }

    // Establish the fixed set of parent summary vectors embedded in the .FLUX
    // file, together with the sampling throttle.
    void initializeFluxSummarySampling_()
    {
        if (this->fluxDumpers_.empty()) {
            return;
        }

        const auto intervalInDays = Parameters::Get<
            Parameters::FluxSummaryMinIntervalBetweenSamples<Scalar>>();

        this->fluxSummaryMinInterval_ =
            unit::convert::from(static_cast<double>(intervalInDays), unit::day);

        const auto boundaryIntervalInDays = Parameters::Get<
            Parameters::FluxBoundaryMinIntervalBetweenSamples<Scalar>>();

        this->fluxBoundaryMinInterval_ =
            unit::convert::from(static_cast<double>(boundaryIntervalInDays), unit::day);

        this->fluxBoundaryReportStepsOnly_ =
            Parameters::Get<Parameters::FluxBoundaryReportStepsOnly>();

        this->fluxBoundaryWindowStart_ = 0.0;
        this->fluxBoundaryHasRecord_ = false;

        for (auto& dumper : this->fluxDumpers_) {
            dumper.setBoundaryMinSampleInterval(this->fluxBoundaryReportStepsOnly_
                                                ? 0.0
                                                : this->fluxBoundaryMinInterval_);
        }

        if (this->fluxBoundaryReportStepsOnly_) {
            OpmLog::note("DUMPFLUX will write sector boundary data once per report step");
        }
        else {
            OpmLog::note(fmt::format("DUMPFLUX will write sector boundary data every time "
                                     "step, at most every {} day(s)",
                                     boundaryIntervalInDays));
        }

        this->fluxSummaryKeyList_ = this->fluxSummaryKeys_();

        this->fluxSummaryKeyTypes_.clear();
        this->fluxSummaryKeyTypes_.reserve(this->fluxSummaryKeyList_.size());
        for (const auto& key : this->fluxSummaryKeyList_) {
            const auto colon = key.find(':');
            const auto keyword = (colon == std::string::npos)
                ? key
                : key.substr(0, colon);

            this->fluxSummaryKeyTypes_.push_back(parseKeywordType(keyword));
        }

        this->fluxSummaryRateAccum_.assign(this->fluxSummaryKeyList_.size(), 0.0);
        this->fluxSummaryAccumDt_ = 0.0;
        this->fluxSummaryLastSampleTime_ = 0.0;
        this->fluxSummaryHasSample_ = false;

        for (auto& dumper : this->fluxDumpers_) {
            dumper.setSummaryKeys(this->fluxSummaryKeyList_);
            dumper.setSummaryMinSampleInterval(this->fluxSummaryMinInterval_);
        }

        OpmLog::note(fmt::format("DUMPFLUX will embed {} summary vector(s), "
                                 "sampled at most every {} day(s)",
                                 this->fluxSummaryKeyList_.size(),
                                 intervalInDays));
    }

    // The FLORES buffers owned by the output module are moved into the restart
    // solution by assignToSolution(). DUMPFLUX therefore has to sample the
    // Cartesian boundary-face values before that happens.
    void captureFluxDumperFlores_()
    {
        this->fluxCapturedFaceRates_.clear();

        if (this->fluxDumpers_.empty()) {
            return;
        }

        const auto& flows = this->outputModule_->getFlows();
        if (!flows.hasFlores()) {
            return;
        }

        this->fluxCapturedFaceRates_.reserve(this->fluxDumpers_.size());
        for (const auto& dumper : this->fluxDumpers_) {
            // The phase volumetric flux across one boundary face, oriented
            // positive into the sector.
            const auto phaseFlux = [&flows, this](const FluxRegions::BoundaryFace& face,
                                                  const EclIO::FluxFile::Phase phase)
            {
                if (face.isNnc || face.direction == FaceDir::Unknown) {
                    return 0.0;
                }

                const auto& vg = this->simulator_.vanguard();

                // Reported by whoever owns the interior cell, once.
                if (!this->fluxCellOwned_(vg.compressedIndex(face.interiorGlobalCell))) {
                    return 0.0;
                }

                const auto comp = this->fluxComponentIndex_(phase);

                // The FLORES buffers are indexed by this rank's own cell
                // numbering, not by cartesian position, so the cartesian index
                // a face carries has to be translated first. Handing the
                // cartesian one straight over reads a different cell entirely
                // as soon as the grid has an inactive cell before this one.
                const auto floresAt = [&flows, &vg, comp](const int globalCell,
                                                          const FaceDir::DirEnum dir)
                {
                    const auto cell = vg.compressedIndex(globalCell);
                    if (cell < 0) {
                        return 0.0;
                    }

                    return flows.getFloresIfAvailable(static_cast<unsigned>(cell), dir, comp);
                };

                // FLORES is stored per cell for the positive face
                // directions only, and is oriented along the positive
                // axis. The FLUX file convention is positive into the
                // sector.
                //
                // For a boundary face on a positive direction the
                // interior cell holds the value and the orientation has
                // to be flipped. For a boundary face on a negative
                // direction the value lives on the exterior cell's
                // positive face and already points into the sector.
                switch (face.direction) {
                case FaceDir::XPlus:
                case FaceDir::YPlus:
                case FaceDir::ZPlus:
                    return -floresAt(face.interiorGlobalCell, face.direction);

                case FaceDir::XMinus:
                    return floresAt(face.exteriorGlobalCell, FaceDir::XPlus);

                case FaceDir::YMinus:
                    return floresAt(face.exteriorGlobalCell, FaceDir::YPlus);

                case FaceDir::ZMinus:
                    return floresAt(face.exteriorGlobalCell, FaceDir::ZPlus);

                default:
                    return 0.0;
                }
            };

            // Converted to component masses here for the same reason the
            // per-time-step path does it: only this side of the boundary knows
            // the state of the cell an inflowing stream comes from.
            this->fluxCapturedFaceRates_.push_back(
                dumper.makeFaceMajorRates(
                    [&phaseFlux, this](const FluxRegions::BoundaryFace& face,
                                       const EclIO::FluxFile::Phase component)
                    {
                        return this->fluxComponentMass_(
                            face, component,
                            [&phaseFlux, &face](const EclIO::FluxFile::Phase p)
                            {
                                return phaseFlux(face, p);
                            });
                    }));
        }
    }

    //! \brief Record the exterior cell's relative permeability and capillary
    //!        pressure for every boundary face of a Pressure-mode record.
    void collectFluxExteriorRockState_(const FluxDumper& dumper,
                                       typename FluxDumper::ReportStepData& step) const
    {
        const auto& vanguard = this->simulator_.vanguard();
        const auto& faces = dumper.data().boundaryFaces;
        const auto numPhaseSlots = static_cast<std::size_t>(dumper.data().header.numPhases);

        step.relPerm.assign(faces.size() * numPhaseSlots, 0.0);
        step.capPressure.assign(faces.size() * numPhaseSlots, 0.0);

        const auto refPhaseIdx = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
            ? FluidSystem::oilPhaseIdx
            : (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)
               ? FluidSystem::gasPhaseIdx
               : FluidSystem::waterPhaseIdx);

        for (std::size_t f = 0; f < faces.size(); ++f) {
            const auto cell = vanguard.compressedIndex(faces[f].exteriorGlobalCell);

            // Reported by the rank that owns the cell and nobody else, so that
            // the caller's sum across ranks picks each face up once.
            if (!this->fluxCellOwned_(cell)) {
                continue;
            }

            const auto* intQuants =
                this->simulator_.model().cachedIntensiveQuantities(cell, /*timeIdx=*/0);
            if (intQuants == nullptr) {
                continue;
            }

            const auto& fs = intQuants->fluidState();

            // The capillary pressures are read off the converged phase
            // pressures rather than from a fresh saturation-function
            // evaluation, so whatever hysteresis state the exterior cell
            // carries is already reflected in them.
            const auto pRef = getValue(fs.pressure(refPhaseIdx));

            std::size_t slot = 0;
            const auto store = [&](const unsigned phaseIdx)
            {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    return;
                }

                step.relPerm[f * numPhaseSlots + slot] =
                    getValue(intQuants->relativePermeability(phaseIdx));
                step.capPressure[f * numPhaseSlots + slot] =
                    getValue(fs.pressure(phaseIdx)) - pRef;
                ++slot;
            };

            // Same canonical Oil/Water/Gas order as the rate arrays.
            store(FluidSystem::oilPhaseIdx);
            store(FluidSystem::waterPhaseIdx);
            store(FluidSystem::gasPhaseIdx);
        }
    }

    //! \brief Pore-volume weighted sums over the cells OUTSIDE a dumper's
    //!        region, in the layout RateConverter accumulates internally.
    //!
    //! \details A reduced run adds these to its own sums so that the averages
    //!   driving a reservoir-volume target are taken over the whole of the
    //!   original model. Sums rather than averages, so that the reduced run
    //!   still reflects changes made inside its own region.
    void collectFluxConverterExternal_(const FluxDumper& dumper,
                                       typename FluxDumper::ReportStepData& step) const
    {
        const auto& vanguard = this->simulator_.vanguard();
        const auto& model = this->simulator_.model();

        std::vector<char> inRegion(model.numGridDof(), 0);
        for (const auto globalCell : dumper.data().localToGlobal) {
            const auto cell = vanguard.compressedIndex(globalCell);
            if (cell >= 0 && static_cast<std::size_t>(cell) < inRegion.size()) {
                inRegion[cell] = 1;
            }
        }

        std::array<double, 8> hpv{};
        std::array<double, 8> pv{};

        const auto pressurePhaseIdx = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
            ? FluidSystem::oilPhaseIdx
            : (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)
               ? FluidSystem::gasPhaseIdx
               : FluidSystem::waterPhaseIdx);

        for (std::size_t cell = 0; cell < inRegion.size(); ++cell) {
            if (inRegion[cell] != 0) {
                continue;
            }

            // Counted by the rank that owns the cell and nobody else. An
            // overlap copy would otherwise add the same pore volume again on
            // every rank that borders it, and these are sums, not averages.
            if (!this->fluxCellOwned_(static_cast<int>(cell))) {
                continue;
            }

            const auto* intQuants = model.cachedIntensiveQuantities(static_cast<unsigned>(cell), 0);
            if (intQuants == nullptr) {
                continue;
            }

            const auto& fs = intQuants->fluidState();
            const double cellPv = model.dofTotalVolume(cell) * getValue(intQuants->porosity());
            if (!(cellPv > 0.0)) {
                continue;
            }

            double hydrocarbon = 1.0;
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                hydrocarbon -= getValue(fs.saturation(FluidSystem::waterPhaseIdx));
            }

            const auto accumulate = [&](std::array<double, 8>& out, const double weight)
            {
                if (!(weight > 0.0)) {
                    return;
                }

                out[0] += getValue(fs.pressure(pressurePhaseIdx)) * weight;
                out[1] += getValue(fs.temperature(pressurePhaseIdx)) * weight;
                out[6] += weight;
                if constexpr (fluxBlackOil_) {
                    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
                        && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))
                    {
                        out[2] += getValue(fs.Rs()) * weight;
                        out[3] += getValue(fs.Rv()) * weight;
                    }
                    if (FluidSystem::enableDissolvedGasInWater()) {
                        out[4] += getValue(fs.Rsw()) * weight;
                    }
                    if (FluidSystem::enableVaporizedWater()) {
                        out[5] += getValue(fs.Rvw()) * weight;
                    }
                    out[7] += getValue(fs.saltConcentration()) * weight;
                }
            };

            accumulate(hpv, cellPv * hydrocarbon);
            accumulate(pv, cellPv);
        }

        step.externalRegionSums.assign(hpv.begin(), hpv.end());
        step.externalRegionSums.insert(step.externalRegionSums.end(), pv.begin(), pv.end());
    }

    // The boundary transmissibilities are not available while the writer is
    // being constructed, so record them the first time a report step is dumped.
    // The consumer uses these to reproduce the parent's inter-cell
    // transmissibility instead of the (larger) default outer-boundary value.
    void assignFluxDumperTransmissibilities_()
    {
        if (this->fluxTransmissibilitiesAssigned_ || this->fluxDumpers_.empty()) {
            return;
        }

        this->fluxTransmissibilitiesAssigned_ = true;

        const auto& vanguard = this->simulator_.vanguard();
        const auto& problem = this->simulator_.problem();

        for (auto& dumper : this->fluxDumpers_) {
            const auto& faces = dumper.regionBoundaryFaces();

            std::vector<double> trans(faces.size(), 0.0);
            std::vector<int> pvtRegion(faces.size(), 0);
            std::vector<double> exteriorDepth(faces.size(), 0.0);
            std::vector<int> exteriorEquilRegion(faces.size(), 0);
            std::vector<double> thermalTrans(faces.size(), 0.0);
            for (std::size_t i = 0; i < faces.size(); ++i) {
                const auto& face = faces[i];
                const auto interior = vanguard.compressedIndex(face.interiorGlobalCell);
                const auto exterior = vanguard.compressedIndex(face.exteriorGlobalCell);
                if (interior < 0 || exterior < 0) {
                    continue;
                }

                // The face belongs to whoever owns its interior cell, so that
                // exactly one rank contributes it to the sum below. The
                // exterior cell need only be reachable, which it is: it borders
                // an owned cell and so sits in this rank's overlap.
                if (!this->fluxCellOwned_(interior)) {
                    continue;
                }

                // Cartesian adjacency does not imply a connection: faults with
                // throw and pinched-out cells leave neighbouring active cells
                // with no shared face, and hence no transmissibility. Such a
                // face carries no flow.
                trans[i] = problem.transmissibilityOrZero(static_cast<unsigned>(interior),
                                                          static_cast<unsigned>(exterior));

                pvtRegion[i] = problem.pvtRegionIndex(static_cast<unsigned>(exterior));

                // A reduced run imposes the exterior cell's temperature at the
                // face, so it needs both half-cells' conduction, not its own half.
                const auto& eclTrans = problem.eclTransmissibilities();
                const double hIn = eclTrans.thermalHalfTransOrZero(static_cast<unsigned>(interior),
                                                                   static_cast<unsigned>(exterior));
                const double hOut = eclTrans.thermalHalfTransOrZero(static_cast<unsigned>(exterior),
                                                                    static_cast<unsigned>(interior));
                if (hIn > 0.0 && hOut > 0.0) {
                    thermalTrans[i] = hIn*hOut / (hIn + hOut);
                }

                // How deep the exterior cell's centre is, which is where the
                // pressures written per report step are taken. A reduced run
                // imposes them at the face, so it needs this to carry them
                // there; in a dipping layer the two are metres apart.
                exteriorDepth[i] = problem.dofCenterDepth(static_cast<unsigned>(exterior));

                // Which equilibration region the exterior cell sits in. A
                // reduced run needs it to find the threshold pressure across
                // the face: its own EQLNUM stops at the boundary, so it has no
                // way of knowing what lies beyond.
                exteriorEquilRegion[i] = static_cast<int>
                    (problem.thresholdPressure()
                     .equilRegionIndex(static_cast<unsigned>(exterior)));
            }

            this->fluxReduceSum_(trans);
            this->fluxReduceSum_(pvtRegion);
            this->fluxReduceSum_(exteriorDepth);
            this->fluxReduceSum_(exteriorEquilRegion);
            this->fluxReduceSum_(thermalTrans);

            // A face nothing flows through says nothing about the cell beyond
            // it, and on more than one rank that cell may not even be at hand:
            // with no connection it is in nobody's overlap, and the sums above
            // leave its entries at zero where a serial run would have filled
            // them in. Record it as unknown whatever the partition, so the file
            // is the same on any number of ranks.
            for (std::size_t i = 0; i < faces.size(); ++i) {
                if (!(trans[i] > 0.0)) {
                    pvtRegion[i] = 0;
                    exteriorDepth[i] = std::numeric_limits<double>::quiet_NaN();
                    exteriorEquilRegion[i] = -1;
                }
            }

            dumper.setBoundaryTransmissibilities(trans);
            dumper.setBoundaryExteriorPvtRegions(pvtRegion);
            dumper.setBoundaryExteriorDepths(exteriorDepth);
            dumper.setBoundaryThermalTransmissibilities(thermalTrans);
            dumper.setBoundaryExteriorEquilRegions(exteriorEquilRegion);

            // The whole region-pair table, not just the entries this boundary
            // happens to touch. A reduced run has interior region boundaries
            // too, and would get those as wrong as the ones at its edge. Every
            // rank holds the same table: it is reduced across them as it is
            // built.
            const auto thpres = problem.thresholdPressure().getRestartVector();
            dumper.setThresholdPressure(std::vector<double>(thpres.begin(), thpres.end()));
        }
    }

    void updateFluxDumpers_(const int reportStepNum, const bool isSubStep)
    {
        if (this->fluxDumpers_.empty()) {
            return;
        }

        // Boundary flow generally varies within a report step, so a record is
        // written per time step by default. The throttle keeps that from
        // producing an unreasonable number of records when the time steps are
        // short; report step boundaries are always written, so the
        // report-step-only behaviour is a strict subset of the default.
        const auto endTime = static_cast<double>(simulator_.time())
            + static_cast<double>(simulator_.timeStepSize());
        const auto elapsed = endTime - this->fluxBoundaryWindowStart_;

        const bool emit = this->fluxBoundaryReportStepsOnly_
            ? !isSubStep
            : (!isSubStep
               || !this->fluxBoundaryHasRecord_
               || (elapsed >= this->fluxBoundaryMinInterval_));

        if (!emit) {
            return;
        }

        this->assignFluxDumperTransmissibilities_();

        const auto& flows = this->outputModule_->getFlows();
        const auto& floresn = this->collectOnIORank_.isParallel()
            ? this->collectOnIORank_.globalFloresn()
            : flows.getFloresn();
        const auto simStep = simulator_.timeStepIndex();

        // The record covers everything since the previous record, which is not
        // the same as the last time step once several steps are accumulated.
        const auto startTime = this->fluxBoundaryWindowStart_;
        const auto stepLength = endTime - this->fluxBoundaryWindowStart_;

        const auto pressurePhaseIdx = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
            ? FluidSystem::oilPhaseIdx
            : (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)
               ? FluidSystem::gasPhaseIdx
               : FluidSystem::waterPhaseIdx);
        const auto& vanguard = this->simulator_.vanguard();

        if (!flows.anyFlores() && !this->fluxMissingFloresReported_) {
            const auto anyFluxMode =
                std::any_of(this->fluxDumpers_.begin(), this->fluxDumpers_.end(),
                            [](const auto& dumper)
                            {
                                return (static_cast<int>(dumper.data().header.mode)
                                        & static_cast<int>(EclIO::FluxFile::Mode::Flux)) != 0;
                            });

            if (anyFluxMode) {
                this->fluxMissingFloresReported_ = true;
                OpmLog::warning("DUMPFLUX is writing FLUX-mode boundary rates, but no FLORES "
                                "data is available. All boundary rates will be zero and a "
                                "USEFLUX run will behave as a closed region. A DUMPFLUX run "
                                "normally turns the FLORES computation on by itself, so this "
                                "means the run was not recognised as a producer.");
            }
        }

        const bool haveCapturedRates =
            this->fluxCapturedFaceRates_.size() == this->fluxDumpers_.size();

        const bool haveAggregatedRates =
            !this->fluxRateTimeWeights_.empty()
            && this->fluxMassSnapshots_.size() == this->fluxDumpers_.size();

        for (std::size_t dumperIdx = 0; dumperIdx < this->fluxDumpers_.size(); ++dumperIdx) {
            auto& dumper = this->fluxDumpers_[dumperIdx];

            auto step = dumper.makeZeroFluxStep(reportStepNum,
                                                simStep,
                                                startTime,
                                                stepLength);

            if (haveAggregatedRates) {
                // Time-weighted average over the time steps that make up this
                // record's window, so that the value a consumer holds across
                // the window is the mass actually transferred over it,
                // independent of how the window was subdivided.
                if (this->fluxMassUsable_
                    && (dumperIdx < this->fluxMassSnapshots_.size())
                    && (this->fluxMassSnapshots_[dumperIdx].size()
                        == this->fluxRateTimeWeights_.size()))
                {
                    step.massRates =
                        FluxDumper::aggregateRates(EclIO::FluxFile::Sampling::Averaged,
                                                   this->fluxMassSnapshots_[dumperIdx],
                                                   this->fluxRateTimeWeights_);
                }
            }
            else if (haveCapturedRates) {
                // Cartesian contributions were sampled before the FLORES buffers
                // were moved into the restart solution.
                step.massRates = this->fluxCapturedFaceRates_[dumperIdx];
            }

            if (haveAggregatedRates || haveCapturedRates) {
                // NNC contributions are not affected by that move and are
                // resolved here so the parallel gather has completed.
                const auto& nncIndex = this->fluxNncPairToIndex_();
                const auto nncRates = dumper.makeFaceMajorRates(
                    [&floresn, &nncIndex, this](const FluxRegions::BoundaryFace& face,
                                                const EclIO::FluxFile::Phase phase)
                    {
                        if (!face.isNnc) {
                            return 0.0;
                        }

                        const auto comp = this->fluxComponentIndex_(phase);
                        const auto key = normalizedNncPair_(face.interiorGlobalCell,
                                                            face.exteriorGlobalCell);
                        const auto it = nncIndex.find(key);
                        if (it == nncIndex.end()) {
                            return 0.0;
                        }

                        const auto nncIdx = static_cast<std::size_t>(it->second);
                        if (nncIdx >= floresn[comp].values.size()) {
                            return 0.0;
                        }

                        // For NNC values we assume the stored direction is cell1->cell2
                        // of the normalized pair and flip sign if the interior is cell1.
                        const auto nncFlux = floresn[comp].values[nncIdx];
                        return (face.interiorGlobalCell == key.second)
                            ? nncFlux
                            : -nncFlux;
                    });

                // NNC faces carry a phase volumetric flux like any other, so it
                // goes through the same upwind and Rs/Rv treatment before being
                // added to the masses.
                if (this->fluxMassUsable_
                    && (step.massRates.size() == nncRates.size()))
                {
                    static constexpr std::array phaseOrder {
                        EclIO::FluxFile::Phase::Oil,
                        EclIO::FluxFile::Phase::Water,
                        EclIO::FluxFile::Phase::Gas,
                    };

                    const auto& header = dumper.data().header;
                    const auto& faces = dumper.regionBoundaryFaces();

                    std::size_t slot = 0;
                    for (const auto& face : faces) {
                        const auto base = slot;

                        // Pick this face's phase rate out of the face-major
                        // vector that was just built.
                        const auto rateOf = [&header, &nncRates, base]
                            (const EclIO::FluxFile::Phase wanted)
                        {
                            std::size_t s = base;
                            for (const auto phase : phaseOrder) {
                                if (!header.hasPhase(phase)) {
                                    continue;
                                }
                                if (phase == wanted) {
                                    return nncRates[s];
                                }
                                ++s;
                            }

                            return 0.0;
                        };

                        for (const auto phase : phaseOrder) {
                            if (!header.hasPhase(phase)) {
                                continue;
                            }
                            if (face.isNnc) {
                                step.massRates[slot] +=
                                    this->fluxComponentMass_(face, phase, rateOf);
                            }
                            ++slot;
                        }
                    }
                }
            }

            if ((static_cast<int>(dumper.data().header.mode) & static_cast<int>(EclIO::FluxFile::Mode::Pressure)) != 0) {
                step.pressures.clear();
                step.swat.clear();
                step.sgas.clear();
                step.rs.clear();
                step.rv.clear();
                step.temperature.clear();

                const auto& faces = dumper.data().boundaryFaces;
                step.pressures.reserve(faces.size());
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    step.swat.reserve(faces.size());
                }
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    step.sgas.reserve(faces.size());
                }
                step.rs.reserve(faces.size());
                step.rv.reserve(faces.size());

                for (const auto& face : faces) {
                    const auto compressedExterior = vanguard.compressedIndex(face.exteriorGlobalCell);

                    // Only the rank that owns the exterior cell reports it, so
                    // that the sum below picks the value up exactly once.
                    // Everyone else leaves zeroes in its place.
                    if (!this->fluxCellOwned_(compressedExterior)) {
                        step.pressures.push_back(0.0);
                        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                            step.swat.push_back(0.0);
                        }
                        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                            step.sgas.push_back(0.0);
                        }
                        step.rs.push_back(0.0);
                        step.rv.push_back(0.0);
                        continue;
                    }

                    const auto& intQuants = *simulator_.model().cachedIntensiveQuantities(compressedExterior, /*timeIdx=*/0);
                    const auto& fs = intQuants.fluidState();
                    step.pressures.push_back(getValue(fs.pressure(pressurePhaseIdx)));
                    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                        step.swat.push_back(getValue(fs.saturation(FluidSystem::waterPhaseIdx)));
                    }
                    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                        step.sgas.push_back(getValue(fs.saturation(FluidSystem::gasPhaseIdx)));
                    }
                    if constexpr (fluxBlackOil_) {
                        step.rs.push_back(getValue(fs.Rs()));
                        step.rv.push_back(getValue(fs.Rv()));
                    }
                    else {
                        step.rs.push_back(0.0);
                        step.rv.push_back(0.0);
                    }
                }

                // Saturation-function state of the exterior cell. The consumer
                // has no cell there, so without these it would have to evaluate
                // its own rock curves at the exterior saturations, using the
                // wrong SATNUM region, the wrong scaled end points and none of
                // this run's hysteresis history.
                this->collectFluxExteriorRockState_(dumper, step);

                // Each rank has filled only the faces whose exterior cell it
                // owns, so collect them. Collective, and reached by every rank
                // because the dumpers exist everywhere.
                this->fluxReduceSum_(step.pressures);
                this->fluxReduceSum_(step.swat);
                this->fluxReduceSum_(step.sgas);
                this->fluxReduceSum_(step.rs);
                this->fluxReduceSum_(step.rv);
                this->fluxReduceSum_(step.relPerm);
                this->fluxReduceSum_(step.capPressure);
            }

            // Written in both modes: a reduced run needs these whatever kind of
            // boundary it uses, because the averages they rebuild drive the
            // wells, not the boundary.
            this->collectFluxConverterExternal_(dumper, step);
            this->fluxReduceSum_(step.externalRegionSums);
            this->fluxReduceSum_(step.massRates);

            // Only the writer keeps the record. The others took part purely to
            // report on the cells they own.
            if (this->collectOnIORank_.isIORank()) {
                dumper.appendReportStep(step);
            }
        }

        // Start a fresh accumulation window for the next record.
        for (auto& snapshots : this->fluxMassSnapshots_) {
            snapshots.clear();
        }
        this->fluxRateTimeWeights_.clear();
        this->fluxBoundaryWindowStart_ = endTime;
        this->fluxBoundaryHasRecord_ = true;

        // Each write appends one block rather than rewriting the file, so this
        // costs the same whether it happens often or rarely. Kept at report
        // step boundaries so a consumer never sees a half-written window.
        if (!isSubStep) {
            this->flushFluxDumpers_();
        }
    }

    void flushFluxDumpers_()
    {
        // One writer only. Every rank holds a dumper so that it can report on
        // its own cells, but the records were gathered onto the IO rank and it
        // is the only one with a file to add them to. Letting the others in
        // here has them all appending to the same path.
        if (!this->collectOnIORank_.isIORank()) {
            return;
        }

        for (std::size_t i = 0; i < this->fluxDumpers_.size(); ++i) {
            this->fluxDumpers_[i].flush(this->fluxOutputPaths_[i], /*formatted=*/false);
        }
    }

    std::vector<std::string> fluxSummaryKeys_() const
    {
        const auto& schedule = this->simulator_.vanguard().schedule();

        // Region, segment and the like cannot be enumerated from the schedule,
        // so hand over the keys this run is configured to produce. That list
        // already covers the vectors an ACTIONX or UDQ needs: opm-common
        // registers those alongside whatever the SUMMARY section asked for,
        // which is why the parent can evaluate a condition on a region its own
        // SUMMARY section never mentions.
        const auto& summaryConfig = this->simulator_.vanguard().summaryConfig();

        auto availableKeys = std::vector<std::string>{};
        availableKeys.reserve(summaryConfig.size());
        for (const auto& node : summaryConfig) {
            availableKeys.push_back(node.uniqueNodeKey());
        }

        return Opm::fluxSummaryKeys(schedule,
                                    FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx),
                                    FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx),
                                    FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx),
                                    availableKeys);
    }

    std::vector<double> fluxSummaryValues_(const std::vector<std::string>& keys) const
    {
        std::vector<double> values;
        values.reserve(keys.size());

        const auto& summaryState = this->summaryState();
        for (const auto& key : keys) {
            values.push_back(summaryState.get(key, 0.0));
        }

        return values;
    }

    // Take one snapshot of the parent summary vectors for the DUMPFLUX file.
    //
    // Called once per time step, after Summary::eval() has populated the
    // SummaryState. Rate-type quantities are accumulated over every step and
    // emitted as a time-average over the interval since the previous retained
    // sample, so that a USEFLUX run holding the value piecewise-constant
    // across that interval reproduces this run's production over it exactly.
    void sampleFluxSummary_(const bool isSubStep, const Scalar sampleTime)
    {
        if (this->fluxDumpers_.empty() || !this->collectOnIORank_.isIORank()) {
            return;
        }

        const auto& keys = this->fluxSummaryKeyList_;
        if (keys.empty()) {
            return;
        }

        const auto values = this->fluxSummaryValues_(keys);
        const auto dt = static_cast<double>(this->simulator_.timeStepSize());

        for (std::size_t i = 0; i < keys.size(); ++i) {
            if (this->fluxSummaryKeyTypes_[i] == SummaryConfigNode::Type::Rate) {
                this->fluxSummaryRateAccum_[i] += values[i] * dt;
            }
        }
        this->fluxSummaryAccumDt_ += dt;

        const auto elapsed = static_cast<double>(sampleTime) - this->fluxSummaryLastSampleTime_;
        const bool retain = !this->fluxSummaryHasSample_
            || !isSubStep
            || (elapsed >= this->fluxSummaryMinInterval_);

        if (!retain) {
            return;
        }

        auto sample = values;
        if (this->fluxSummaryAccumDt_ > 0.0) {
            for (std::size_t i = 0; i < keys.size(); ++i) {
                if (this->fluxSummaryKeyTypes_[i] == SummaryConfigNode::Type::Rate) {
                    sample[i] = this->fluxSummaryRateAccum_[i] / this->fluxSummaryAccumDt_;
                }
            }
        }

        for (auto& dumper : this->fluxDumpers_) {
            dumper.appendSummarySample(static_cast<double>(sampleTime), sample);
        }

        std::fill(this->fluxSummaryRateAccum_.begin(), this->fluxSummaryRateAccum_.end(), 0.0);
        this->fluxSummaryAccumDt_ = 0.0;
        this->fluxSummaryLastSampleTime_ = static_cast<double>(sampleTime);
        this->fluxSummaryHasSample_ = true;
    }

    // Fill in the rates of wells that lie outside the USEFLUX region from the
    // parent run, so that summary aggregation covers the whole field.
    //
    // Three conventions matter here: data::Rates holds production as NEGATIVE
    // values, it is expressed in SI units while summary vectors are in the
    // deck's output units, and Summary::eval() skips any well whose dynamic
    // status is SHUT - which is exactly what a well with no active connection
    // would otherwise be.
    void injectParentOnlyWellData_(data::Wells& wellData, const Scalar time)
    {
        const auto* parent = this->simulator_.problem().fluxParentSummary();
        if (parent == nullptr) {
            return;
        }

        const auto& vanguard = this->simulator_.vanguard();
        const auto& schedule = vanguard.schedule();
        const auto& units = vanguard.eclState().getUnits();
        const auto stepIdx =
            static_cast<std::size_t>(std::max(this->simulator_.episodeIndex(), 0));

        // The rates go into this step's cumulatives, so they have to be the
        // parent's average over the interval since the summary was last
        // evaluated, however many of the parent's own samples that spans.
        const auto start = (this->fluxParentOutputTime_ < 0.0)
            ? static_cast<double>(time)
            : this->fluxParentOutputTime_;
        this->fluxParentOutputTime_ = static_cast<double>(time);
        this->fluxParentOutputStart_ = start;

        for (const auto& wellName : schedule.wellNames(stepIdx)) {
            const auto& well = schedule.getWell(wellName, stepIdx);

            const auto hasLocalConnection =
                std::any_of(well.getConnections().begin(), well.getConnections().end(),
                            [&vanguard](const auto& conn)
                            {
                                return vanguard.compressedIndex(conn.global_index()) >= 0;
                            });

            if (hasLocalConnection) {
                continue;
            }

            const auto isProducer = well.isProducer();
            const auto sign = isProducer ? -1.0 : 1.0;

            auto& target = wellData[wellName];
            target.dynamicStatus = Well::Status::OPEN;

            const auto assign = [&](const data::Rates::opt opt,
                                    const UnitSystem::measure measure,
                                    const std::string& prodKeyword,
                                    const std::string& injKeyword)
            {
                const auto key = (isProducer ? prodKeyword : injKeyword) + ':' + wellName;
                const auto value = parent->valueOver(key, start, static_cast<double>(time));

                if (std::isfinite(value)) {
                    target.rates.set(opt, sign * units.to_si(measure, value));
                }
            };

            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                assign(data::Rates::opt::oil, UnitSystem::measure::liquid_surface_rate,
                       "WOPR", "WOIR");
            }
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                assign(data::Rates::opt::wat, UnitSystem::measure::liquid_surface_rate,
                       "WWPR", "WWIR");
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                assign(data::Rates::opt::gas, UnitSystem::measure::gas_surface_rate,
                       "WGPR", "WGIR");
            }

            // The parent reports one reservoir volume per well, so it is booked
            // where the well model books it for a well outside the sector: a
            // producer's against oil, or the first active phase without oil,
            // and an injector's against the phase it injects.
            const auto reservoirOpt = [&well]()
            {
                if (well.isProducer()) {
                    return FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
                        ? data::Rates::opt::reservoir_oil
                        : (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
                           ? data::Rates::opt::reservoir_water
                           : data::Rates::opt::reservoir_gas);
                }
                switch (well.injectorType()) {
                case InjectorType::GAS: return data::Rates::opt::reservoir_gas;
                case InjectorType::OIL: return data::Rates::opt::reservoir_oil;
                default:                return data::Rates::opt::reservoir_water;
                }
            }();
            assign(reservoirOpt, UnitSystem::measure::rate, "WVPR", "WVIR");

            // Pressures are not rates and are taken at the end of the interval.
            // A UDQ may read them, and without this Summary::eval() reports
            // zero for a well it has no solution for, over the parent's value.
            const auto pressure = [&](const std::string& keyword, double& target_value)
            {
                const auto value = parent->valueAt(keyword + ':' + wellName, static_cast<double>(time));
                if (std::isfinite(value)) {
                    target_value = units.to_si(UnitSystem::measure::pressure, value);
                }
            };
            pressure("WBHP", target.bhp);
            pressure("WTHP", target.thp);
        }
    }

    //! \brief Replace with the parent's what a sector cannot compute, before
    //!        the UDQs are evaluated.
    //!
    //! \details Called between Summary::eval() and the UDQ evaluation. The
    //!   parent's values were seeded into the SummaryState at the start of
    //!   the step, but Summary::eval() then overwrites every vector the run
    //!   reports with what it computes itself, and for anything that lies
    //!   outside the sector that is zero, or covers the sector's cells only: a
    //!   well it has no solution for, a segment of one, a region it holds part
    //!   of, a block it does not have, a field pressure. A UDQ over any of
    //!   those, and every ACTIONX and group target that reads the UDQ, would
    //!   then go its own way. So these are put back to the parent's values,
    //!   which also makes them what the run reports.
    //!
    //!   What the sector does compute correctly is left alone: every quantity
    //!   of a well it simulates, the rates and totals of field and group from
    //!   which the parent-only wells' rates are already summed, cumulative
    //!   quantities, which accumulate from those rates, and the UDQs, which
    //!   are evaluated next.
    void applyParentOnlySummaryValues_(SummaryState& summaryState, const Scalar time)
    {
        const auto* parent = this->simulator_.problem().fluxParentSummary();
        if (parent == nullptr) {
            return;
        }

        // A UDQ over a region the sector does not hold, such as RGPR_EQL 1
        // with all of EQLNUM 1 outside it, is refused when its region set
        // has no such region. The regions the parent reported are the ones
        // its UDQs read, so those are made known.
        if (!this->udqRegionStatistics_.has_value()) {
            auto statistics = this->eclState_.fipRegionStatistics();
            for (const auto& key : parent->keys()) {
                if (const auto region = ParentSummary::regionOf(key); region.has_value()) {
                    statistics.raiseMaximumRegionID(region->first, region->second);
                }
            }
            this->udqRegionStatistics_ = std::move(statistics);
        }

        const auto& vanguard = this->simulator_.vanguard();
        const auto& schedule = vanguard.schedule();
        const auto& grid = vanguard.eclState().getInputGrid();
        const auto stepIdx =
            static_cast<std::size_t>(std::max(this->simulator_.episodeIndex(), 0));

        // In a USEFLUX run the schedule was built on the sector's cells, so a
        // well outside it has no connections, whether the deck completed it
        // outside the sector or the well was added from the FLUX file.
        std::unordered_set<std::string> absentWells;
        for (const auto& wellName : schedule.wellNames(stepIdx)) {
            if (schedule.getWell(wellName, stepIdx).getConnections().empty()) {
                absentWells.insert(wellName);
            }
        }

        // Given to the parent-only wells by injectParentOnlyWellData_(), or
        // summed over them by Summary::eval().
        static const std::unordered_set<std::string> fromWellRates {
            "OPR", "WPR", "GPR", "VPR", "LPR", "OIR", "WIR", "GIR", "VIR", "LIR",
            "GOR", "WCT", "GLR", "OGR", "WGR", "BHP", "THP",
        };

        const auto isUdq = [](const std::string& keyword)
        {
            return (keyword.size() > 1) && (keyword[1] == 'U');
        };

        const auto start = this->fluxParentOutputStart_;
        const auto end = static_cast<double>(time);

        for (const auto& key : parent->keys()) {
            const auto colon = key.find(':');
            const auto keyword = key.substr(0, colon);
            const auto type = ParentSummary::keyType(key);

            if (isUdq(keyword) || (type == SummaryConfigNode::Type::Total)) {
                continue;
            }

            const auto quantity = keyword.substr(1);

            bool fromParent = false;
            switch (EclIO::SummaryNode::category_from_keyword(keyword)) {
            case EclIO::SummaryNode::Category::Well:
            case EclIO::SummaryNode::Category::Connection:
            case EclIO::SummaryNode::Category::Completion:
            case EclIO::SummaryNode::Category::Segment: {
                if (colon == std::string::npos) {
                    break;
                }
                const auto second = key.find(':', colon + 1);
                const auto wellName = key.substr(colon + 1, (second == std::string::npos)
                                                 ? std::string::npos : second - colon - 1);
                const auto isWellVector = EclIO::SummaryNode::category_from_keyword(keyword)
                    == EclIO::SummaryNode::Category::Well;
                fromParent = (absentWells.count(wellName) > 0)
                    && !(isWellVector && (fromWellRates.count(quantity) > 0));
                break;
            }

            case EclIO::SummaryNode::Category::Field:
            case EclIO::SummaryNode::Category::Group:
                fromParent = fromWellRates.count(quantity) == 0;
                break;

            case EclIO::SummaryNode::Category::Region:
                fromParent = true;
                break;

            case EclIO::SummaryNode::Category::Block: {
                int i = 0, j = 0, k = 0;
                if ((colon != std::string::npos)
                    && (std::sscanf(key.c_str() + colon + 1, "%d,%d,%d", &i, &j, &k) == 3))
                {
                    fromParent = !grid.cellActive(i - 1, j - 1, k - 1);
                }
                break;
            }

            default:
                break;
            }

            if (!fromParent) {
                continue;
            }

            const auto value = (type == SummaryConfigNode::Type::Rate)
                ? parent->valueOver(key, start, end)
                : parent->valueAt(key, end);

            if (std::isfinite(value)) {
                ParentSummary::assign(summaryState, key, value);
            }
        }
    }

    Simulator& simulator_;
    std::unique_ptr<OutputModule> outputModule_;
    Scalar restartTimeStepSize_;
    int rank_ ;
    Inplace inplace_;
    std::vector<FluxDumper> fluxDumpers_;
    std::vector<std::string> fluxOutputPaths_;
    //! Position of each NNC in the writer's NNC list, keyed by its normalised
    //! Cartesian pair. Built on first use, since the list is only complete
    //! once the transmissibilities have been exported. See fluxNncPairToIndex_().
    mutable std::map<std::pair<int, int>, int> fluxNncIndex_;
    std::vector<std::vector<double>> fluxCapturedFaceRates_;
    std::vector<std::vector<std::vector<double>>> fluxMassSnapshots_;
    std::vector<double> fluxRateTimeWeights_;
    std::vector<std::string> fluxSummaryKeyList_;
    std::vector<SummaryConfigNode::Type> fluxSummaryKeyTypes_;
    std::vector<double> fluxSummaryRateAccum_;
    double fluxSummaryAccumDt_ = 0.0;
    //! When the parent-only wells were last given rates, or negative before
    //! the first time. See injectParentOnlyWellData_().
    double fluxParentOutputTime_ = -1.0;
    //! Start of the interval the parent-only wells were last given rates for.
    double fluxParentOutputStart_ = 0.0;
    double fluxSummaryLastSampleTime_ = 0.0;
    double fluxSummaryMinInterval_ = 0.0;
    bool fluxSummaryHasSample_ = false;
    double fluxBoundaryMinInterval_ = 0.0;
    double fluxBoundaryWindowStart_ = 0.0;
    bool fluxBoundaryHasRecord_ = false;
    bool fluxBoundaryReportStepsOnly_ = false;
    bool fluxMassUsable_ = false;
    bool fluxMissingFloresReported_ = false;

    //! \brief Set once the boundary component masses have been reported as
    //!        disagreeing with the assembled residual.
    bool fluxMassMismatchReported_ = false;
    bool fluxTransmissibilitiesAssigned_ = false;

    //! \brief Whether this rank owns each local cell, i.e. holds it as an
    //!        interior cell rather than as an overlap copy.
    std::vector<char> fluxOwnedCell_;
};

} // namespace Opm

#endif // OPM_ECL_WRITER_HPP
