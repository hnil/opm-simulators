/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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


#ifndef OPM_BLACKOILWELLMODEL_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_HEADER_INCLUDED

#include <ebos/eclproblem.hh>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <cassert>
#include <tuple>

#include <opm/parser/eclipse/EclipseState/Runspec.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTestState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GuideRate.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GConSale.hpp>

#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/simulators/flow/countGlobalCells.hpp>
#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/WellInterface.hpp>
#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/simulators/wells/MultisegmentWell.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/timestepping/gatherConvergenceReport.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmatrix.hh>

#include <opm/material/densead/Math.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>

BEGIN_PROPERTIES

NEW_PROP_TAG(EnableTerminalOutput);

END_PROPERTIES

namespace Opm {

        /// Class for handling the blackoil well model.
        template<typename TypeTag>
        class BlackoilWellModel : public Opm::BaseAuxiliaryModule<TypeTag>
        {
        public:
            // ---------      Types      ---------
            typedef WellStateFullyImplicitBlackoil WellState;
            typedef BlackoilModelParametersEbos<TypeTag> ModelParameters;

            typedef typename GET_PROP_TYPE(TypeTag, Grid)                Grid;
            typedef typename GET_PROP_TYPE(TypeTag, FluidSystem)         FluidSystem;
            typedef typename GET_PROP_TYPE(TypeTag, ElementContext)      ElementContext;
            typedef typename GET_PROP_TYPE(TypeTag, Indices)             Indices;
            typedef typename GET_PROP_TYPE(TypeTag, Simulator)           Simulator;
            typedef typename GET_PROP_TYPE(TypeTag, Scalar)              Scalar;
            typedef typename GET_PROP_TYPE(TypeTag, RateVector)          RateVector;
            typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector)      GlobalEqVector;
            typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;

            typedef typename Opm::BaseAuxiliaryModule<TypeTag>::NeighborSet NeighborSet;

            static const int numEq = Indices::numEq;
            static const int solventSaturationIdx = Indices::solventSaturationIdx;

            // TODO: where we should put these types, WellInterface or Well Model?
            // or there is some other strategy, like TypeTag
            typedef Dune::FieldVector<Scalar, numEq    > VectorBlockType;
            typedef Dune::BlockVector<VectorBlockType> BVector;

            typedef Dune::FieldMatrix<Scalar, numEq, numEq > MatrixBlockType;

            typedef Opm::BlackOilPolymerModule<TypeTag> PolymerModule;

            // For the conversion between the surface volume rate and resrevoir voidage rate
            using RateConverterType = RateConverter::
                SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;

            BlackoilWellModel(Simulator& ebosSimulator);

            void init();

            /////////////
            // <eWoms auxiliary module stuff>
            /////////////
            unsigned numDofs() const
            // No extra dofs are inserted for wells. (we use a Schur complement.)
            { return 0; }

            void addNeighbors(std::vector<NeighborSet>& neighbors) const;

            void applyInitial()
            {}

            void linearize(SparseMatrixAdapter& jacobian, GlobalEqVector& res);

            void postSolve(GlobalEqVector& deltaX)
            {
                recoverWellSolutionAndUpdateWellState(deltaX);
            }

            /////////////
            // </ eWoms auxiliary module stuff>
            /////////////

            template <class Restarter>
            void deserialize(Restarter& /* res */)
            {
                // TODO (?)
            }

            /*!
             * \brief This method writes the complete state of the well
             *        to the harddisk.
             */
            template <class Restarter>
            void serialize(Restarter& /* res*/)
            {
                // TODO (?)
            }

            void beginEpisode()
            {
                beginReportStep(ebosSimulator_.episodeIndex());
            }

            void beginTimeStep();

            void beginIteration()
            {
                assemble(ebosSimulator_.model().newtonMethod().numIterations(),
                         ebosSimulator_.timeStepSize());
            }

            void solveWells(double dt){
                last_report_ = SimulatorReportSingle();
                if ( ! wellsActive() ) {
                    return;
                }
                bool converged;
                int exception_thrown = 0;
                Opm::DeferredLogger local_deferredLogger;
                try{
                    updatePerforationIntensiveQuantities();
                    std::vector< Scalar > B_avg(numComponents(), Scalar() );
                    computeAverageFormationFactor(B_avg);
                    updatePerforationIntensiveQuantities();
                    last_report_ = solveWellEq(B_avg, dt, local_deferredLogger);
                } catch (std::exception& e) {
                    exception_thrown = 1;
                }
                logAndCheckForExceptionsAndThrow(local_deferredLogger, exception_thrown, "assemble() failed.", terminal_output_);
                
                last_report_.converged = true;
            }
            
            void endIteration()
            { }

            void endTimeStep()
            {
                timeStepSucceeded(ebosSimulator_.time(), ebosSimulator_.timeStepSize());
            }

            void endEpisode()
            {
                endReportStep();
            }

            template <class Context>
            void computeTotalRatesForDof(RateVector& rate,
                                         const Context& context,
                                         unsigned spaceIdx,
                                         unsigned timeIdx) const;


            using WellInterfacePtr = std::shared_ptr<WellInterface<TypeTag> >;
            WellInterfacePtr well(const std::string& wellName) const;

            void initFromRestartFile(const RestartValue& restartValues);

            Opm::data::Group groupData(const int reportStepIdx, Opm::Schedule& sched) const
            {
                Opm::data::Group dw;
                for (const std::string gname :  sched.groupNames(reportStepIdx))  {
                    const auto& grup = sched.getGroup(gname, reportStepIdx);
                    const auto& grup_type = grup.getGroupType();
                    Opm::data::currentGroupConstraints cgc;
                    cgc.currentProdConstraint =  Opm::Group::ProductionCMode::NONE;
                    cgc.currentGasInjectionConstraint = Opm::Group::InjectionCMode::NONE;
                    cgc.currentWaterInjectionConstraint = Opm::Group::InjectionCMode::NONE;
                    if (this->well_state_.hasProductionGroupControl(gname)) {
                        cgc.currentProdConstraint = this->well_state_.currentProductionGroupControl(gname);
                    }
                    if ((grup_type == Opm::Group::GroupType::INJECTION) || (grup_type == Opm::Group::GroupType::MIXED))  {
                        if (this->well_state_.hasInjectionGroupControl(Opm::Phase::WATER, gname)) {
                            cgc.currentWaterInjectionConstraint = this->well_state_.currentInjectionGroupControl(Opm::Phase::WATER, gname);
                        }
                        if (this->well_state_.hasInjectionGroupControl(Opm::Phase::GAS, gname)) {
                            cgc.currentGasInjectionConstraint = this->well_state_.currentInjectionGroupControl(Opm::Phase::GAS, gname);
                        }
                    }
                    dw.emplace(gname, cgc);
                }
                return dw;
            }

            Opm::data::Wells wellData() const
            {
                auto wsrpt = well_state_.report(phase_usage_, Opm::UgGridHelpers::globalCell(grid()));

                for (const auto& well : this->wells_ecl_) {
                    auto xwPos = wsrpt.find(well.name());
                    if (xwPos == wsrpt.end()) { // No well results.  Unexpected.
                        continue;
                    }

                    xwPos->second.current_control.isProducer = well.isProducer();
                }

                return wsrpt;
            }

            // substract Binv(D)rw from r;
            void apply( BVector& r) const;

            // subtract B*inv(D)*C * x from A*x
            void apply(const BVector& x, BVector& Ax) const;

#if HAVE_CUDA
            // accumulate the contributions of all Wells in the WellContributions object
            void getWellContributions(WellContributions& x) const;
#endif

            // apply well model with scaling of alpha
            void applyScaleAdd(const Scalar alpha, const BVector& x, BVector& Ax) const;

            // Check if well equations is converged.
            ConvergenceReport getWellConvergence(const std::vector<Scalar>& B_avg, const bool checkGroupConvergence = false) const;

            // return the internal well state, ignore the passed one.
            // Used by the legacy code to make it compatible with the legacy well models.
            const WellState& wellState(const WellState& well_state OPM_UNUSED) const;

            // return the internal well state
            const WellState& wellState() const;

            // return the internal well state
            const WellState& prevWellState() const;

            const SimulatorReportSingle& lastReport() const;

            void addWellContributions(SparseMatrixAdapter& jacobian) const
            {
                for ( const auto& well: well_container_ ) {
                    well->addWellContributions(jacobian);
                }
            }

            // called at the beginning of a report step
            void beginReportStep(const int time_step);

            /// Return true if any well has a THP constraint.
            bool hasTHPConstraints() const;

            /// Shut down any single well, but only if it is in prediction mode.
            /// Returns true if the well was actually found and shut.
            bool forceShutWellByNameIfPredictionMode(const std::string& wellname, const double simulation_time);

        protected:
            Simulator& ebosSimulator_;

            std::vector< Well > wells_ecl_;
            std::vector< std::vector<PerforationData> > well_perf_data_;
            std::vector<int> first_perf_index_;

            bool wells_active_;

            // a vector of all the wells.
            std::vector<WellInterfacePtr > well_container_;

            // map from logically cartesian cell indices to compressed ones
            std::vector<int> cartesian_to_compressed_;

            std::vector<bool> is_cell_perforated_;

            void initializeWellPerfData();

            // create the well container
            std::vector<WellInterfacePtr > createWellContainer(const int time_step);

            WellInterfacePtr createWellForWellTest(const std::string& well_name, const int report_step, Opm::DeferredLogger& deferred_logger) const;

            WellState well_state_;
            WellState previous_well_state_;
            WellState well_state_nupcol_;

            const ModelParameters param_;
            bool terminal_output_;
            bool has_solvent_;
            bool has_polymer_;
            std::vector<int> pvt_region_idx_;
            PhaseUsage phase_usage_;
            size_t global_nc_;
            // the number of the cells in the local grid
            size_t number_of_cells_;
            double gravity_;
            std::vector<double> depth_;
            bool initial_step_;
            bool report_step_starts_;

            std::unique_ptr<RateConverterType> rateConverter_;
            std::unique_ptr<VFPProperties<VFPInjProperties,VFPProdProperties>> vfp_properties_;

            SimulatorReportSingle last_report_;

            WellTestState wellTestState_;
            std::unique_ptr<GuideRate> guideRate_;

            // used to better efficiency of calcuation
            mutable BVector scaleAddRes_;

            const Grid& grid() const
            { return ebosSimulator_.vanguard().grid(); }

            const EclipseState& eclState() const
            { return ebosSimulator_.vanguard().eclState(); }

            const Schedule& schedule() const
            { return ebosSimulator_.vanguard().schedule(); }

            // compute the well fluxes and assemble them in to the reservoir equations as source terms
            // and in the well equations.
            void assemble(const int iterationIdx,
                          const double dt);

            // called at the end of a time step
            void timeStepSucceeded(const double& simulationTime, const double dt);

            // called at the end of a report step
            void endReportStep();

            // using the solution x to recover the solution xw for wells and applying
            // xw to update Well State
            void recoverWellSolutionAndUpdateWellState(const BVector& x);

            void updateWellControls(Opm::DeferredLogger& deferred_logger, const bool checkGroupControls);

            void updateAndCommunicateGroupData();

            // setting the well_solutions_ based on well_state.
            void updatePrimaryVariables(Opm::DeferredLogger& deferred_logger);

            void setupCartesianToCompressed_(const int* global_cell, int number_of_cells);

            void computeRepRadiusPerfLength(const Grid& grid, Opm::DeferredLogger& deferred_logger);


            void computeAverageFormationFactor(std::vector<Scalar>& B_avg) const;

            // Calculating well potentials for each well
            void computeWellPotentials(std::vector<double>& well_potentials, const int reportStepIdx, Opm::DeferredLogger& deferred_logger);

            const std::vector<double>& wellPerfEfficiencyFactors() const;

            void calculateEfficiencyFactors(const int reportStepIdx);

            // it should be able to go to prepareTimeStep(), however, the updateWellControls() and initPrimaryVariablesEvaluation()
            // makes it a little more difficult. unless we introduce if (iterationIdx != 0) to avoid doing the above functions
            // twice at the beginning of the time step
            /// Calculating the explict quantities used in the well calculation. By explicit, we mean they are cacluated
            /// at the beginning of the time step and no derivatives are included in these quantities
            void calculateExplicitQuantities(Opm::DeferredLogger& deferred_logger) const;

            SimulatorReportSingle solveWellEq(const std::vector<Scalar>& B_avg, const double dt, Opm::DeferredLogger& deferred_logger);

            void initPrimaryVariablesEvaluation() const;

            // The number of components in the model.
            int numComponents() const;

            int numLocalWells() const;

            int numPhases() const;

            void assembleWellEq(const std::vector<Scalar>& B_avg, const double dt, Opm::DeferredLogger& deferred_logger);

            // some preparation work, mostly related to group control and RESV,
            // at the beginning of each time step (Not report step)
            void prepareTimeStep(Opm::DeferredLogger& deferred_logger);

            void extractLegacyCellPvtRegionIndex_();

            void extractLegacyDepth_();

            /// return true if wells are available in the reservoir
            bool wellsActive() const;

            void setWellsActive(const bool wells_active);

            /// return true if wells are available on this process
            bool localWellsActive() const;

            /// upate the wellTestState related to economic limits
            void updateWellTestState(const double& simulationTime, WellTestState& wellTestState) const;

            void updatePerforationIntensiveQuantities();

            void wellTesting(const int timeStepIdx, const double simulationTime, Opm::DeferredLogger& deferred_logger);

            // convert well data from opm-common to well state from opm-core
            void wellsToState( const data::Wells& wells,
                               const PhaseUsage& phases,
                               const bool handle_ms_well,
                               WellStateFullyImplicitBlackoil& state ) const;

            // whether there exists any multisegment well open on this process
            bool anyMSWellOpenLocal() const;

            const Well& getWellEcl(const std::string& well_name) const;

            void updateGroupIndividualControls(Opm::DeferredLogger& deferred_logger, std::set<std::string>& switched_groups);
            void updateGroupIndividualControl(const Group& group, Opm::DeferredLogger& deferred_logger, std::set<std::string>& switched_groups);
            bool checkGroupConstraints(const Group& group, Opm::DeferredLogger& deferred_logger) const;
            Group::ProductionCMode checkGroupProductionConstraints(const Group& group, Opm::DeferredLogger& deferred_logger) const;
            Group::InjectionCMode checkGroupInjectionConstraints(const Group& group, const Phase& phase) const;
            void checkGconsaleLimits(const Group& group, WellState& well_state, Opm::DeferredLogger& deferred_logger ) const;

            void updateGroupHigherControls(Opm::DeferredLogger& deferred_logger, std::set<std::string>& switched_groups);
            void checkGroupHigherConstraints(const Group& group, Opm::DeferredLogger& deferred_logger, std::set<std::string>& switched_groups);

            void actionOnBrokenConstraints(const Group& group, const Group::ExceedAction& exceed_action, const Group::ProductionCMode& newControl, Opm::DeferredLogger& deferred_logger);

            void actionOnBrokenConstraints(const Group& group, const Group::InjectionCMode& newControl, const Phase& topUpPhase, Opm::DeferredLogger& deferred_logger);

            WellInterfacePtr getWell(const std::string& well_name) const;

            void updateWsolvent(const Group& group, const Schedule& schedule, const int reportStepIdx, const WellStateFullyImplicitBlackoil& wellState);

            void setWsolvent(const Group& group, const Schedule& schedule, const int reportStepIdx, double wsolvent);



        };


} // namespace Opm

#include "BlackoilWellModel_impl.hpp"
#endif
