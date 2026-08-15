/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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
#include <opm/simulators/flow/BlackoilModelParameters.hpp>

#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/nonlinear/newtonmethodparams.hpp>

#include <opm/models/utils/parametersystem.hpp>

#include <algorithm>
#include <stdexcept>

namespace Opm {

template<class Scalar>
BlackoilModelParameters<Scalar>::BlackoilModelParameters()
{
    dbhp_max_rel_ = Parameters::Get<Parameters::DbhpMaxRel<Scalar>>();
    dwell_fraction_max_ = Parameters::Get<Parameters::DwellFractionMax<Scalar>>();
    inj_mult_osc_threshold_ = Parameters::Get<Parameters::InjMultOscThreshold<Scalar>>();
    inj_mult_damp_mult_ = Parameters::Get<Parameters::InjMultDampMult<Scalar>>();
    inj_mult_min_damp_factor_ = Parameters::Get<Parameters::InjMultMinDampFactor<Scalar>>();
    max_residual_allowed_ = Parameters::Get<Parameters::MaxResidualAllowed<Scalar>>();
    relaxed_max_pv_fraction_ = Parameters::Get<Parameters::RelaxedMaxPvFraction<Scalar>>();
    cnv_pv_floor_fraction_ = Parameters::Get<Parameters::CnvPoreVolumeFloorFraction<Scalar>>();
    cnv_throughput_weight_ = Parameters::Get<Parameters::CnvThroughputWeight<Scalar>>();
    tolerance_mb_ = Parameters::Get<Parameters::ToleranceMb<Scalar>>();
    tolerance_mb_relaxed_ = std::max(tolerance_mb_, Parameters::Get<Parameters::ToleranceMbRelaxed<Scalar>>());
    tolerance_energy_balance_ = Parameters::Get<Parameters::ToleranceEnergyBalance<Scalar>>();
    tolerance_energy_balance_relaxed_ = std::max(tolerance_energy_balance_, Parameters::Get<Parameters::ToleranceEnergyBalanceRelaxed<Scalar>>());
    tolerance_cnv_ = Parameters::Get<Parameters::ToleranceCnv<Scalar>>();
    tolerance_cnv_relaxed_ = std::max(tolerance_cnv_, Parameters::Get<Parameters::ToleranceCnvRelaxed<Scalar>>());
    tolerance_cnv_energy_ = Parameters::Get<Parameters::ToleranceCnvEnergy<Scalar>>();
    tolerance_cnv_energy_relaxed_ = std::max(tolerance_cnv_energy_, Parameters::Get<Parameters::ToleranceCnvEnergyRelaxed<Scalar>>());
    tolerance_max_dp_ = Parameters::Get<Parameters::ToleranceMaxDp<Scalar>>();
    tolerance_max_ds_ = Parameters::Get<Parameters::ToleranceMaxDs<Scalar>>();
    tolerance_max_drs_ = Parameters::Get<Parameters::ToleranceMaxDrs<Scalar>>();
    tolerance_max_drv_ = Parameters::Get<Parameters::ToleranceMaxDrv<Scalar>>();
    tolerance_wells_ = Parameters::Get<Parameters::ToleranceWells<Scalar>>();
    tolerance_wells_stopped_factor_ = Parameters::Get<Parameters::ToleranceWellsStoppedFactor<Scalar>>();
    tolerance_wells_dynamic_thp_factor_ = Parameters::Get<Parameters::ToleranceWellsDynamicThpFactor<Scalar>>();
    max_welleq_iter_ = Parameters::Get<Parameters::MaxWelleqIter>();
    use_multisegment_well_ = Parameters::Get<Parameters::UseMultisegmentWell>();
    tolerance_pressure_ms_wells_ = Parameters::Get<Parameters::TolerancePressureMsWells<Scalar>>();
    tolerance_well_bhp_eq_ = Parameters::Get<Parameters::ToleranceWellBhpEq<Scalar>>();
    tolerance_well_thp_eq_ = Parameters::Get<Parameters::ToleranceWellThpEq<Scalar>>();
    tolerance_well_grup_eq_ = Parameters::Get<Parameters::ToleranceWellGrupEq<Scalar>>();
    well_thp_switch_band_ = Parameters::Get<Parameters::WellThpSwitchBand<Scalar>>();
    well_vfp_explicit_reset_ = Parameters::Get<Parameters::WellVfpExplicitReset>();
    relaxed_tolerance_flow_well_ = Parameters::Get<Parameters::RelaxedWellFlowTol<Scalar>>();
    relaxed_tolerance_pressure_ms_well_ = Parameters::Get<Parameters::RelaxedPressureTolMsw<Scalar>>();
    max_pressure_change_ms_wells_ = Parameters::Get<Parameters::MaxPressureChangeMsWells<Scalar>>();
    max_inner_iter_ms_wells_ = Parameters::Get<Parameters::MaxInnerIterMsWells>();
    strict_inner_iter_wells_ = Parameters::Get<Parameters::StrictInnerIterWells>();
    strict_outer_iter_wells_ = Parameters::Get<Parameters::StrictOuterIterWells>();
    regularization_factor_wells_ = Parameters::Get<Parameters::RegularizationFactorWells<Scalar>>();
    max_niter_inner_well_iter_ = Parameters::Get<Parameters::MaxNewtonIterationsWithInnerWellIterations>();
    shut_unsolvable_wells_ = Parameters::Get<Parameters::ShutUnsolvableWells>();
    max_inner_iter_wells_ = Parameters::Get<Parameters::MaxInnerIterWells>();
    min_iter_after_switch_wells_ = Parameters::Get<Parameters::MinIterAfterSwitchWells>();
    min_iter_after_switch_ms_wells_ = Parameters::Get<Parameters::MinIterAfterSwitchMsWells>();
    max_well_status_switch_inner_iter_ = Parameters::Get<Parameters::MaxWellStatusSwitchInInnerIterWells>();
    max_well_status_switch_ = Parameters::Get<Parameters::MaxWellStatusSwitchForWells>();
    maxSinglePrecisionTimeStep_ = Parameters::Get<Parameters::MaxSinglePrecisionDays<Scalar>>() * 24 * 60 * 60;
    min_strict_cnv_iter_ = Parameters::Get<Parameters::MinStrictCnvIter>();
    min_strict_mb_iter_ = Parameters::Get<Parameters::MinStrictMbIter>();
    solve_welleq_initially_ = Parameters::Get<Parameters::SolveWelleqInitially>();
    pre_solve_network_ = Parameters::Get<Parameters::PreSolveNetwork>();
    update_equations_scaling_ = Parameters::Get<Parameters::UpdateEquationsScaling>();
    use_update_stabilization_ = Parameters::Get<Parameters::UseUpdateStabilization>();
    matrix_add_well_contributions_ = Parameters::Get<Parameters::MatrixAddWellContributions>();
    check_well_operability_ = Parameters::Get<Parameters::EnableWellOperabilityCheck>();
    check_well_operability_iter_ = Parameters::Get<Parameters::EnableWellOperabilityCheckIter>();
    max_number_of_well_switches_ = Parameters::Get<Parameters::MaximumNumberOfWellSwitches>();
    well_control_flip_flop_hold_ = Parameters::Get<Parameters::WellControlFlipFlopHold>();
    well_inner_control_pin_ = Parameters::Get<Parameters::WellInnerControlPin>();
    well_inner_line_search_ = Parameters::Get<Parameters::WellInnerLineSearch>();
    well_inner_line_search_start_ = Parameters::Get<Parameters::WellInnerLineSearchStart>();
    well_inner_line_search_armijo_ = Parameters::Get<Parameters::WellInnerLineSearchArmijo<Scalar>>();
    well_inner_line_search_min_step_ = Parameters::Get<Parameters::WellInnerLineSearchMinStep<Scalar>>();
    well_inner_damping_ = Parameters::Get<Parameters::WellInnerDamping>();
    well_inner_min_relaxation_ = Parameters::Get<Parameters::WellInnerMinRelaxation<Scalar>>();
    well_inner_stall_exit_ = Parameters::Get<Parameters::WellInnerStallExit>();
    well_inner_cycle_exit_ = Parameters::Get<Parameters::WellInnerCycleExit>();
    well_inner_return_best_ = Parameters::Get<Parameters::WellInnerReturnBest>();
    well_inner_frozen_exit_ = Parameters::Get<Parameters::WellInnerFrozenExit>();
    well_status_hysteresis_ = Parameters::Get<Parameters::WellStatusHysteresis<Scalar>>();
    well_status_keep_on_cap_ = Parameters::Get<Parameters::WellStatusKeepOnCap>();
    well_inner_trace_ = Parameters::Get<Parameters::WellInnerTrace>();
    newton_confirm_iterations_ = Parameters::Get<Parameters::NewtonConfirmIterations>();
    gate_ensemble_subiterate_ = Parameters::Get<Parameters::GateEnsembleSubIterate>();
    max_number_of_group_switches_ = Parameters::Get<Parameters::MaximumNumberOfGroupSwitches>();
    use_average_density_ms_wells_ = Parameters::Get<Parameters::UseAverageDensityMsWells>();
    local_well_solver_control_switching_ = Parameters::Get<Parameters::LocalWellSolveControlSwitching>();
    use_implicit_ipr_ = Parameters::Get<Parameters::UseImplicitIpr>();
    check_group_constraints_inner_well_iterations_ = Parameters::Get<Parameters::CheckGroupConstraintsInnerWellIterations>();
    nonlinear_solver_ = Parameters::Get<Parameters::NonlinearSolver>();
    const auto approach = Parameters::Get<Parameters::LocalSolveApproach>();
    if (approach == "jacobi") {
        local_solve_approach_ = DomainSolveApproach::Jacobi;
    } else if (approach == "gauss-seidel") {
        local_solve_approach_ = DomainSolveApproach::GaussSeidel;
    } else {
        throw std::runtime_error("Invalid domain solver approach '" + approach + "' specified.");
    }

    max_local_solve_iterations_ = Parameters::Get<Parameters::MaxLocalSolveIterations>();
    local_tolerance_scaling_mb_ = Parameters::Get<Parameters::LocalToleranceScalingMb<Scalar>>();
    local_tolerance_scaling_cnv_ = Parameters::Get<Parameters::LocalToleranceScalingCnv<Scalar>>();
    newton_max_iter_ = Parameters::Get<Parameters::NewtonMaxIterations>();
    newton_min_iter_ = Parameters::Get<Parameters::NewtonMinIterations>();
    nldd_num_initial_newton_iter_ = Parameters::Get<Parameters::NlddNumInitialNewtonIter>();
    nldd_relative_mobility_change_tol_ = Parameters::Get<Parameters::NlddRelativeMobilityChangeTol<Scalar>>();
    num_local_domains_ = Parameters::Get<Parameters::NumLocalDomains>();
    local_domains_partition_imbalance_ = std::max(Scalar{1.0}, Parameters::Get<Parameters::LocalDomainsPartitioningImbalance<Scalar>>());
    local_domains_partition_method_ = Parameters::Get<Parameters::LocalDomainsPartitioningMethod>();
    local_domains_partition_well_neighbor_levels_ = Parameters::Get<Parameters::LocalDomainsPartitionWellNeighborLevels>();
    deck_file_name_ = Parameters::Get<Parameters::EclDeckFileName>();
    network_max_strict_outer_iterations_ = Parameters::Get<Parameters::NetworkMaxStrictOuterIterations>();
    network_max_outer_iterations_ = Parameters::Get<Parameters::NetworkMaxOuterIterations>();
    network_max_sub_iterations_ = Parameters::Get<Parameters::NetworkMaxSubIterations>();
    rc_network_loose_coupling_ = Parameters::Get<Parameters::RcNetworkLooseCoupling>();
    network_pressure_update_damping_factor_ = Parameters::Get<Parameters::NetworkPressureUpdateDampingFactor<Scalar>>();
    network_max_pressure_update_in_bars_ = Parameters::Get<Parameters::NetworkMaxPressureUpdateInBars<Scalar>>();
    local_domains_ordering_ = domainOrderingMeasureFromString(Parameters::Get<Parameters::LocalDomainsOrderingMeasure>());
    write_partitions_ = Parameters::Get<Parameters::DebugEmitCellPartition>();

    monitor_params_.enabled_ = Parameters::Get<Parameters::ConvergenceMonitoring>();
    monitor_params_.cutoff_ = Parameters::Get<Parameters::ConvergenceMonitoringCutOff>();
    monitor_params_.decay_factor_ = Parameters::Get<Parameters::ConvergenceMonitoringDecayFactor<Scalar>>();

    nupcol_group_rate_tolerance_ = Parameters::Get<Parameters::NupcolGroupRateTolerance<Scalar>>();
    well_group_constraints_max_iterations_ = Parameters::Get<Parameters::WellGroupConstraintsMaxIterations>();
    group_control_fraction_tolerance_ = Parameters::Get<Parameters::GroupControlFractionTolerance<Scalar>>();
}

template<class Scalar>
void BlackoilModelParameters<Scalar>::registerParameters()
{
    Parameters::Register<Parameters::DbhpMaxRel<Scalar>>
        ("Maximum relative change of the bottom-hole pressure in a single iteration");
    Parameters::Register<Parameters::DwellFractionMax<Scalar>>
        ("Maximum absolute change of a well's volume fraction in a single iteration");
    Parameters::Register<Parameters::InjMultOscThreshold<Scalar>>
        ("Injection multiplier oscillation threshold (used for multiplier dampening)");
    Parameters::Register<Parameters::InjMultDampMult<Scalar>>
        ("Injection multiplier dampening factor (dampening multiplied by this each time oscillation is detected)");
    Parameters::Register<Parameters::InjMultMinDampFactor<Scalar>>
        ("Minimum injection multiplier dampening factor (maximum dampening level)");
    Parameters::Register<Parameters::MaxResidualAllowed<Scalar>>
        ("Absolute maximum tolerated for residuals without cutting the time step size");
    Parameters::Register<Parameters::CnvThroughputWeight<Scalar>>
        ("Normalise the CNV residual of a cell by max(PV, w*B*dt*q) instead of the pore "
         "volume alone, with w the value given here and q the sum of the component's flux "
         "magnitudes through the cell (faces, boundary faces, sources and wells). A cell "
         "whose storage is negligible against what flows through it is then judged on "
         "relative flux imbalance, |R| <= tol*w*q, which stays finite as PV -> 0 and is "
         "also the correct measure in the incompressible/quasi-steady limit. Since the "
         "sum of face-flux magnitudes counts both what enters and what leaves, w = 0.5 "
         "corresponds to one pass through the cell. 0 (default) disables, giving the "
         "storage-only criterion.");
    Parameters::Register<Parameters::CnvPoreVolumeFloorFraction<Scalar>>
        ("Floor the per-cell pore volume used in the CNV normalisation at this fraction "
         "of the average cell pore volume. The CNV criterion divides the residual by the "
         "cell pore volume, so cells with vanishing PV (inner radial cells, salt-clogged "
         "cells, degenerate geometry) are held to a divergently strict standard; the "
         "floor bounds the criterion for such cells. 0 (default) disables.");
    Parameters::Register<Parameters::RelaxedMaxPvFraction<Scalar>>
        ("The fraction of the pore volume of the reservoir "
         "where the volumetric error (CNV) may be violated "
         "during strict Newton iterations.");
    Parameters::Register<Parameters::ToleranceMb<Scalar>>
        ("Tolerated mass balance error relative to total mass present");
    Parameters::Register<Parameters::ToleranceMbRelaxed<Scalar>>
        ("Relaxed tolerated mass balance error that applies for iterations "
         "after the iterations with the strict tolerance");
    Parameters::Register<Parameters::ToleranceEnergyBalance<Scalar>>
        ("Tolerated energy balance error relative to (scaled) total energy present");
    Parameters::Register<Parameters::ToleranceEnergyBalanceRelaxed<Scalar>>
        ("Relaxed tolerated energy balance error that applies for iterations "
         "after the iterations with the strict tolerance");
    Parameters::Register<Parameters::ToleranceCnv<Scalar>>
        ("Local convergence tolerance (Maximum of local saturation errors)");
    Parameters::Register<Parameters::ToleranceCnvRelaxed<Scalar>>
        ("Relaxed local convergence tolerance that applies for iterations "
         "after the iterations with the strict tolerance");
    Parameters::Register<Parameters::ToleranceCnvEnergy<Scalar>>
        ("Local energy convergence tolerance (Maximum of local energy errors)");
    Parameters::Register<Parameters::ToleranceCnvEnergyRelaxed<Scalar>>
        ("Relaxed local energy convergence tolerance that applies for iterations "
         "after the iterations with the strict tolerance");
    Parameters::Register<Parameters::ToleranceMaxDp<Scalar>>
        ("Tolerance for max pressure change during a Newton iteration. "
         "A value greater than 0.0 allows for convergence regardless "
         "of residual tolerances. Use with care!");
    Parameters::Register<Parameters::ToleranceMaxDs<Scalar>>
        ("Tolerance for max saturation change during a Newton iteration. "
         "A value greater than 0.0 allows for convergence regardless "
         "of residual tolerances. Use with care!");
    Parameters::Register<Parameters::ToleranceMaxDrs<Scalar>>
        ("Tolerance for max RS change during a Newton iteration. "
         "A value greater than 0.0 allows for convergence regardless "
         "of residual tolerances. Use with care!");
    Parameters::Register<Parameters::ToleranceMaxDrv<Scalar>>
        ("Tolerance for max RV change during a Newton iteration. "
         "A value greater than 0.0 allows for convergence regardless "
         "of residual tolerances. Use with care!");
    Parameters::Register<Parameters::ToleranceWells<Scalar>>
        ("Well convergence tolerance");
    Parameters::Register<Parameters::ToleranceWellsStoppedFactor<Scalar>>
        ("Multiplier applied to the well convergence tolerance for stopped wells "
         "and wells under a zero rate target");
    Parameters::Register<Parameters::ToleranceWellsDynamicThpFactor<Scalar>>
        ("Multiplier applied to the well convergence tolerance for wells on a "
         "dynamic THP limit, to help network convergence");
    Parameters::Register<Parameters::MaxWelleqIter>
        ("Maximum number of iterations to determine solution the well equations");
    Parameters::Register<Parameters::UseMultisegmentWell>
        ("Use the well model for multi-segment wells instead of the "
         "one for single-segment wells");
    Parameters::Register<Parameters::TolerancePressureMsWells<Scalar>>
        ("Tolerance for the pressure equations for multi-segment wells");
    Parameters::Register<Parameters::ToleranceWellBhpEq<Scalar>>
        ("Tolerance for the BHP control equation residual of standard wells [Pa]");
    Parameters::Register<Parameters::ToleranceWellThpEq<Scalar>>
        ("Tolerance for the THP control equation residual of standard wells [Pa]");
    Parameters::Register<Parameters::ToleranceWellGrupEq<Scalar>>
        ("Tolerance for the group-control equation residual of standard wells [m3/s]");
    Parameters::Register<Parameters::WellThpSwitchBand<Scalar>>
        ("Dead band [Pa] on the switch-to-THP decision of individual wells: a well on "
         "rate/BHP control only switches to THP once its THP limit exceeds the current "
         "THP by this much. Both sides of that comparison are approximations (the "
         "current THP is itself a VFP evaluation), and the bare '>' lets a well that "
         "just left THP control re-enter it on tolerance-level noise. 0 (default) keeps "
         "the bare comparison");
    Parameters::Register<Parameters::WellVfpExplicitReset>
        ("Reset the well's explicit-VFP-lookup fallback whenever its operability status "
         "is reset (start of every timestep attempt), instead of letting it latch for "
         "the rest of the timestep once any solve engaged it. Off by default");
    Parameters::Register<Parameters::RelaxedWellFlowTol<Scalar>>
        ("Relaxed tolerance for the well flow residual");
    Parameters::Register<Parameters::RelaxedPressureTolMsw<Scalar>>
        ("Relaxed tolerance for the MSW pressure solution");
    Parameters::Register<Parameters::MaxPressureChangeMsWells<Scalar>>
        ("Maximum relative pressure change for a single iteration "
         "of the multi-segment well model");
    Parameters::Register<Parameters::MaxInnerIterMsWells>
        ("Maximum number of inner iterations for multi-segment wells");
    Parameters::Register<Parameters::StrictInnerIterWells>
        ("Number of inner well iterations with strict tolerance");
    Parameters::Register<Parameters::StrictOuterIterWells>
        ("Number of newton iterations for which wells are checked with strict tolerance");
    Parameters::Register<Parameters::MaxNewtonIterationsWithInnerWellIterations>
        ("Maximum newton iterations with inner well iterations");
    Parameters::Register<Parameters::ShutUnsolvableWells>
        ("Shut unsolvable wells");
    Parameters::Register<Parameters::MaxInnerIterWells>
        ("Maximum number of inner iterations for standard wells");
    Parameters::Register<Parameters::MinIterAfterSwitchWells>
        ("Inner iterations forced between two control switches of a standard well. "
         "Values <= 1 let the inner solve cycle between controls.");
    Parameters::Register<Parameters::MinIterAfterSwitchMsWells>
        ("Inner iterations forced between two control switches of a multi-segment well. "
         "Values <= 1 let the inner solve cycle between controls.");
    Parameters::Register<Parameters::MaxWellStatusSwitchInInnerIterWells>
        ("Maximum number of status switching (stop<->open) for a well during inner iterations.");
    Parameters::Register<Parameters::MaxWellStatusSwitchForWells>
        ("Maximum number of status switching (stop<->open) for a well during a time-step.");
    Parameters::Register<Parameters::AlternativeWellRateInit>
        ("Use alternative well rate initialization procedure");
    Parameters::Register<Parameters::RegularizationFactorWells<Scalar>>
        ("Regularization factor for wells");
    Parameters::Register<Parameters::MaxSinglePrecisionDays<Scalar>>
        ("Maximum time step size where single precision floating point "
         "arithmetic can be used solving for the linear systems of equations");
    Parameters::Register<Parameters::MinStrictCnvIter>
        ("Minimum number of Newton iterations before relaxed tolerances "
         "can be used for the CNV convergence criterion. "
         "Default -1 means that the relaxed tolerance is used when maximum "
         "number of Newton iterations are reached.");
    Parameters::Register<Parameters::MinStrictMbIter>
        ("Minimum number of Newton iterations before relaxed tolerances "
         "can be used for the MB convergence criterion. "
         "Default -1 means that the relaxed tolerance is used when maximum "
         "number of Newton iterations are reached.");
    Parameters::Register<Parameters::SolveWelleqInitially>
        ("Fully solve the well equations before each iteration of the reservoir model");
    Parameters::Register<Parameters::PreSolveNetwork>
        ("Pre solve and iterate the network model at start-up");
    Parameters::Register<Parameters::UpdateEquationsScaling>
        ("Update scaling factors for mass balance equations during the run");
    Parameters::Register<Parameters::UseUpdateStabilization>
        ("Try to detect and correct oscillations or stagnation during the Newton method");
    Parameters::Register<Parameters::MatrixAddWellContributions>
        ("Explicitly specify the influences of wells between cells in "
         "the Jacobian and preconditioner matrices");
    Parameters::Register<Parameters::EnableWellOperabilityCheck>
        ("Enable the well operability checking");
    Parameters::Register<Parameters::EnableWellOperabilityCheckIter>
        ("Enable the well operability checking during iterations");
    Parameters::Register<Parameters::MaximumNumberOfWellSwitches>
        ("Maximum number of times a well can switch to the same control");
    Parameters::Register<Parameters::WellControlFlipFlopHold>
        ("Hold a well's control mode for the remainder of the timestep attempt "
         "after this many immediate control reversals (switching back to the "
         "control it just left). Switches to group control are never refused. "
         "Values below 2 hold before the constraint scan has settled and are "
         "not recommended. 0 (default) disables the hold");
    Parameters::Register<Parameters::WellInnerControlPin>
        ("Pin a well's control mode in the inner well solve after this many "
         "consecutive outer Newton iterations whose inner solve switched the "
         "control away and back (a net-zero flip pair). The well solution then "
         "stops changing control between outer iterations; the outer control "
         "update keeps its authority. 0 (default) disables the pin");
    Parameters::Register<Parameters::WellInnerLineSearch>
        ("Backtracking line search on the standard-well inner Newton step: maximum "
         "number of step halvings per inner iteration, accepting the first step that "
         "reduces the scaled well residual norm. The trial assembly is reused as the "
         "next iteration's assembly, so an accepted full step costs nothing extra. "
         "0 (default) disables the line search");
    Parameters::Register<Parameters::WellInnerLineSearchStart>
        ("Inner iteration from which the line search is engaged. Wells that converge "
         "before it never pay a trial assembly. 0 (default) engages it immediately");
    Parameters::Register<Parameters::WellInnerLineSearchArmijo<Scalar>>
        ("Sufficient-decrease constant of the inner line search: a trial step is "
         "accepted when the merit drops below (1 - c*step) times its current value");
    Parameters::Register<Parameters::WellInnerLineSearchMinStep<Scalar>>
        ("Shortest step fraction the inner line search will try before accepting the "
         "undamped step and moving on");
    Parameters::Register<Parameters::WellInnerDamping>
        ("Damp the standard-well inner Newton step when its residual measure "
         "oscillates, the way multisegment wells already do. Off by default");
    Parameters::Register<Parameters::WellInnerMinRelaxation<Scalar>>
        ("Floor for the oscillation damping factor of the standard-well inner solve");
    Parameters::Register<Parameters::WellInnerStallExit>
        ("Abandon the standard-well inner solve after this many consecutive iterations "
         "in which the line search found no descent direction, instead of spending the "
         "remaining iteration budget. Requires the line search. 0 (default) disables");
    Parameters::Register<Parameters::WellInnerCycleExit>
        ("Leave the standard-well inner solve once it returns for the Nth time to a "
         "(control mode, scaled residual) state it has already been in. The switching "
         "logic is then cycling and no amount of further iteration will converge it. "
         "0 (default) disables the check");
    Parameters::Register<Parameters::WellInnerFrozenExit>
        ("Leave the standard-well inner solve after this many consecutive iterations "
         "in which the iterate did not change at all. 0 disables.");
    Parameters::Register<Parameters::WellStatusHysteresis<Scalar>>
        ("Relative dead band between the well stop and reopen decisions. The reopen "
         "test must clear the operable bhp limit by this fraction, and the well is "
         "reopened inside the operable region rather than on its boundary. 0 disables.");
    Parameters::Register<Parameters::WellStatusKeepOnCap>
        ("When the inner status-switch cap binds, keep the current well status instead "
         "of reverting to the status the inner solve started with.");
    Parameters::Register<Parameters::WellInnerReturnBest>
        ("When the standard-well inner solve ends without converging, return the "
         "lowest-residual iterate it visited rather than wherever it stopped. Off by default");
    Parameters::Register<Parameters::WellInnerTrace>
        ("Log one line per standard-well inner iteration for the named well "
         "(\"*\" for all wells): merit, bhp, rate, control mode, damping. Diagnostic only");
    Parameters::Register<Parameters::GateEnsembleSubIterate>
        ("When an iteration is physically converged but blocked only by the well-group or "
         "network gate, run up to this many sweeps of {solve all wells, refresh group "
         "data} at frozen reservoir state to let the control ensemble settle, instead of "
         "buying full Newton iterations for it. 0 (default) disables.");
    Parameters::Register<Parameters::NewtonConfirmIterations>
        ("Take this many additional full Newton iterations after the convergence "
         "criterion is first met. The accepted state then carries a residual well "
         "below the tolerance instead of just under it, which reduces how far two "
         "nearly-identical runs can drift apart at each accepted step. Costs the "
         "extra iterations on every substep. 0 (default) disables");
    Parameters::Register<Parameters::MaximumNumberOfGroupSwitches>
        ("Maximum number of times a group can switch to the same control");
    Parameters::Register<Parameters::UseAverageDensityMsWells>
        ("Approximate segment densitities by averaging over segment and its outlet");
    Parameters::Register<Parameters::LocalWellSolveControlSwitching>
        ("Allow control switching during local well solutions");
    Parameters::Register<Parameters::UseImplicitIpr>
        ("Compute implict IPR for stability checks and stable solution search");
    Parameters::Register<Parameters::CheckGroupConstraintsInnerWellIterations>
        ("Allow checking of group constraints during inner well iterations");
    Parameters::Register<Parameters::NetworkMaxStrictOuterIterations>
        ("Maximum outer iterations in network solver before relaxing tolerance");
    Parameters::Register<Parameters::NetworkMaxOuterIterations>
        ("Maximum outer number of iterations in the network solver before giving up");
    Parameters::Register<Parameters::NetworkMaxSubIterations>
        ("Maximum number of sub-iterations to update network pressures (within a single well/group control update)");
    Parameters::Register<Parameters::RcNetworkLooseCoupling>
        ("Reservoir coupling: exchange master/slave network pressures and rates only once per master outer "
         "network iteration (loose coupling) instead of the default once per inner sub-iteration (tight coupling)");
    Parameters::Register<Parameters::NetworkPressureUpdateDampingFactor<Scalar>>
        ("Damping factor in the inner network pressure update iterations");
    Parameters::Register<Parameters::NetworkMaxPressureUpdateInBars<Scalar>>
        ("Maximum pressure update in the inner network pressure update iterations");
    Parameters::Register<Parameters::NonlinearSolver>
        ("Choose nonlinear solver. Valid choices are newton or nldd.");
    Parameters::Register<Parameters::LocalSolveApproach>
        ("Choose local solve approach. Valid choices are jacobi and gauss-seidel");
    Parameters::SetDefault<Parameters::NewtonMaxIterations>(20);
    Parameters::Register<Parameters::NewtonMinIterations>
        ("The minimum number of Newton iterations per time step");
    Parameters::Register<Parameters::MaxLocalSolveIterations>
        ("Max iterations for local solves with NLDD nonlinear solver.");
    Parameters::Register<Parameters::LocalToleranceScalingMb<Scalar>>
        ("Set lower than 1.0 to use stricter convergence tolerance for local solves.");
    Parameters::Register<Parameters::LocalToleranceScalingCnv<Scalar>>
        ("Set lower than 1.0 to use stricter convergence tolerance for local solves.");
    Parameters::Register<Parameters::NlddNumInitialNewtonIter>
        ("Number of initial global Newton iterations when running the NLDD nonlinear solver.");
    Parameters::Register<Parameters::NlddRelativeMobilityChangeTol<Scalar>>
        ("Threshold for single cell relative mobility change in the NLDD solver");
    Parameters::Register<Parameters::NumLocalDomains>
        ("Number of local domains for NLDD nonlinear solver.");
    Parameters::Register<Parameters::LocalDomainsPartitioningImbalance<Scalar>>
        ("Subdomain partitioning imbalance tolerance. 1.03 is 3 percent imbalance.");
    Parameters::Register<Parameters::LocalDomainsPartitioningMethod>
        ("Subdomain partitioning method. Allowed values are "
         "'zoltan', "
         "'simple', "
         "and the name of a partition file ending with '.partition'.");
    Parameters::Register<Parameters::LocalDomainsPartitionWellNeighborLevels>
        ("Number of neighbor levels around wells to include in the same domain during NLDD partitioning");
    Parameters::Register<Parameters::LocalDomainsOrderingMeasure>
        ("Subdomain ordering measure. Allowed values are "
         "'maxpressure', "
         "'averagepressure' "
         "and  'residual'.");
    Parameters::Register<Parameters::DebugEmitCellPartition>
        ("Whether or not to emit cell partitions as a debugging aid.");

    Parameters::Register<Parameters::ConvergenceMonitoring>
        ("Enable convergence monitoring");
    Parameters::Register<Parameters::ConvergenceMonitoringCutOff>
        ("Cut off limit for convergence monitoring");
    Parameters::Register<Parameters::ConvergenceMonitoringDecayFactor<Scalar>>
        ("Decay factor for convergence monitoring");

    Parameters::Register<Parameters::NupcolGroupRateTolerance<Scalar>>
        ("Tolerance for acceptable changes in VREP/RAIN group rates");

    Parameters::Hide<Parameters::DebugEmitCellPartition>();

    Parameters::Register<Parameters::WellGroupConstraintsMaxIterations>
        ("Maximum number of iterations in the well/group switching algorithm");

    Parameters::Register<Parameters::GroupControlFractionTolerance<Scalar>>
        ("Tolerance for minimal allowed fraction for controlled phase before employing a fallback control");
    // if openMP is available, use two threads per mpi rank by default
#if _OPENMP
    Parameters::SetDefault<Parameters::ThreadsPerProcess>(2);
#endif
}

template struct BlackoilModelParameters<double>;

#if FLOW_INSTANTIATE_FLOAT
template struct BlackoilModelParameters<float>;
#endif

} // namespace Opm
