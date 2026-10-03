opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-test.sh "")

if (BUILD_FLOW)
  opm_add_executable(
    TARGET
      flux_smoke_shape_check
    SOURCES
      tests/flux/flux_smoke_shape_check.cpp
    LIBRARIES
      opmcommon
  )

  opm_add_executable(
    TARGET
      flux_smoke_mutate
    SOURCES
      tests/flux/flux_smoke_mutate.cpp
    LIBRARIES
      opmcommon
  )

  opm_add_executable(
    TARGET
      flux_smoke_compare
    SOURCES
      tests/flux/flux_smoke_compare.cpp
    LIBRARIES
      opmcommon
  )

  # The summary reader is built and shipped by opm-common, not by
  # opm-simulators.
  find_program(OPM_SUMMARY_BIN
    NAMES summary
    HINTS ${PROJECT_SOURCE_DIR}/../opm-common/build/bin
          ${CMAKE_INSTALL_PREFIX}/bin
    PATH_SUFFIXES bin)

  # compareECL is built and shipped by opm-common, not by opm-simulators.
  find_program(OPM_COMPARE_ECL_BIN
    NAMES compareECL
    HINTS ${PROJECT_SOURCE_DIR}/../opm-common/build/bin
          ${CMAKE_INSTALL_PREFIX}/bin
    PATH_SUFFIXES bin)

  # ------------------------------------------------------------------------
  # DUMPFLUX: a run writing the boundary of its FLUXNUM regions.

  set(_flux_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_dump_smoke)
  opm_add_test(flux_dump_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_flux_result_path}
      -f FLUX_DUMP_SMOKE
      -p "${PROJECT_BINARY_DIR}/bin/flux_smoke_shape_check FLUX_DUMP_SMOKE.FLUX 10 1 10 10 1 5 FLUX"
  )

  set(_flux_multi_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_dump_multi_smoke)
  opm_add_test(flux_dump_multi_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_flux_multi_result_path}
      -f FLUX_DUMP_MULTI_SMOKE
      -p "${PROJECT_BINARY_DIR}/bin/flux_smoke_shape_check FLUX_DUMP_MULTI_SMOKE.FLUX0001 10 1 10 10 1 5 FLUX FLUX_DUMP_MULTI_SMOKE.FLUX0002 10 1 10 1 1 10 FLUX"
  )

  set(_flux_pressure_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_dump_pressure_smoke)
  opm_add_test(flux_dump_pressure_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_flux_pressure_result_path}
      -f FLUX_DUMP_PRESSURE_SMOKE
      -p "${PROJECT_BINARY_DIR}/bin/flux_smoke_shape_check FLUX_DUMP_PRESSURE_SMOKE.FLUX 10 1 10 1 1 2 PRESSURE"
  )

  set(_flux_both_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_dump_both_smoke)
  opm_add_test(flux_dump_both_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_flux_both_result_path}
      -f FLUX_DUMP_BOTH_SMOKE
      -p "${PROJECT_BINARY_DIR}/bin/flux_smoke_shape_check FLUX_DUMP_BOTH_SMOKE.FLUX 10 1 10 10 1 5 BOTH"
  )

  set(_flux_invalid_type_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_dump_invalid_type_smoke)
  opm_add_test(flux_dump_invalid_type_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_flux_invalid_type_result_path}
      -f FLUX_DUMP_INVALID_TYPE_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-flux-invalid-type-smoke.sh ${PROJECT_BINARY_DIR}/bin/flux_smoke_shape_check"
  )

  # FLUXNUM without DUMPFLUX is an ordinary run.
  set(_flux_fluxnum_only_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_fluxnum_only_smoke)
  opm_add_test(flux_fluxnum_only_smoke
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_flux_fluxnum_only_result_path}
      -f FLUX_FLUXNUM_ONLY_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-flux-fluxnum-only-smoke.sh"
  )

  set(_flux_mixed_case_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_dump_mixed_case_smoke)
  opm_add_test(flux_dump_mixed_case_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_flux_mixed_case_result_path}
      -f FLUX_DUMP_MIXED_CASE_SMOKE
      -p "${PROJECT_BINARY_DIR}/bin/flux_smoke_shape_check FLUX_DUMP_MIXED_CASE_SMOKE.FLUX 10 1 10 10 1 5 BOTH"
  )

  # A FLUX file must come out the same however many ranks wrote it. Covered
  # for both payloads: pressure mode carries the exterior state and the
  # region sums, flux mode the component masses.
  foreach(_flux_par_case FLUX_DUMP_PRESSURE_SMOKE FLUX_DUMP_FLORES_SMOKE)
    string(TOLOWER ${_flux_par_case} _flux_par_lower)
    set(_flux_par_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_parallel_${_flux_par_lower})
    opm_add_test(flux_parallel_${_flux_par_lower}
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/flux
        -r ${_flux_par_result_path}
        -f ${_flux_par_case}
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-flux-parallel-consistency.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_BINARY_DIR}/bin/flux_smoke_compare ${PROJECT_SOURCE_DIR}/tests/flux/${_flux_par_case}.DATA 3"
    )
  endforeach()

  # ------------------------------------------------------------------------
  # make_flux: the same boundary data, formed from a finished run's output.

  set(_make_flux_pressure_equivalence_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+make_flux_pressure_equivalence_smoke)
  opm_add_test(make_flux_pressure_equivalence_smoke
    DEPENDS
      flux_smoke_compare
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_make_flux_pressure_equivalence_result_path}
      -f FLUX_DUMP_PRESSURE_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-make-flux-equivalence-smoke.sh ${PROJECT_BINARY_DIR}/bin/make_flux ${PROJECT_BINARY_DIR}/bin/flux_smoke_compare ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_PRESSURE_SMOKE.DATA FLUX_DUMP_PRESSURE_SMOKE pressure 10"
  )

  set(_make_flux_pressure_equivalence_mpi_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+make_flux_pressure_equivalence_mpi_smoke)
  opm_add_test(make_flux_pressure_equivalence_mpi_smoke
    DEPENDS
      flux_smoke_compare
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_make_flux_pressure_equivalence_mpi_result_path}
      -f FLUX_DUMP_PRESSURE_SMOKE
      -n 2
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-make-flux-equivalence-smoke.sh ${PROJECT_BINARY_DIR}/bin/make_flux ${PROJECT_BINARY_DIR}/bin/flux_smoke_compare ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_PRESSURE_SMOKE.DATA FLUX_DUMP_PRESSURE_SMOKE pressure 10"
  )

  set(_make_flux_udq_keys_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+make_flux_udq_keys_smoke)
  opm_add_test(make_flux_udq_keys_smoke
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_make_flux_udq_keys_result_path}
      -f FLUX_DUMP_UDQKEYS_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-make-flux-udq-keys.sh ${PROJECT_BINARY_DIR}/bin/make_flux ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_UDQKEYS_SMOKE.DATA ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_SMOKE_ONE_SECTOR.grdecl"
  )

  set(_make_flux_sparse_restart_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+make_flux_sparse_restart_smoke)
  opm_add_test(make_flux_sparse_restart_smoke
      EXE_TARGET flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_make_flux_sparse_restart_result_path}
      -f FLUX_DUMP_SPARSE_RESTART_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-make-flux-sparse-restart.sh ${PROJECT_BINARY_DIR}/bin/make_flux ${PROJECT_BINARY_DIR}/bin/inspect_flux ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_SPARSE_RESTART_SMOKE.DATA ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_SMOKE_ONE_SECTOR.grdecl"
  )

  set(_make_flux_fluxall_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+make_flux_fluxall_smoke)
  opm_add_test(make_flux_fluxall_smoke
      EXE_TARGET flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_make_flux_fluxall_result_path}
      -f FLUX_DUMP_FLUXALL_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-make-flux-fluxall.sh ${PROJECT_BINARY_DIR}/bin/make_flux ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_FLUXALL_SMOKE.DATA ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_SMOKE_ONE_SECTOR.grdecl"
  )

  set(_make_flux_fluxnum_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+make_flux_fluxnum_smoke)
  opm_add_test(make_flux_fluxnum_smoke
    DEPENDS
      flux_smoke_compare
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_make_flux_fluxnum_result_path}
      -f FLUX_DUMP_PRESSURE_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-make-flux-fluxnum-smoke.sh ${PROJECT_BINARY_DIR}/bin/make_flux ${PROJECT_BINARY_DIR}/bin/flux_smoke_compare ${PROJECT_BINARY_DIR}/bin/flux_smoke_shape_check ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_PRESSURE_SMOKE.DATA ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_SMOKE_ONE_SECTOR.grdecl ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_SMOKE_TWO_SECTORS.grdecl"
  )

  set(_flux_minpv_boundary_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_minpv_boundary_smoke)
  opm_add_test(flux_minpv_boundary_smoke
    DEPENDS
      flux_smoke_compare
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_flux_minpv_boundary_result_path}
      -f FLUX_DUMP_MINPV_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-flux-minpv-boundary-smoke.sh ${PROJECT_BINARY_DIR}/bin/make_flux ${PROJECT_BINARY_DIR}/bin/inspect_flux ${PROJECT_BINARY_DIR}/bin/flux_smoke_compare ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_MINPV_SMOKE.DATA ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_SMOKE_ONE_SECTOR.grdecl"
  )

  set(_make_flux_flux_equivalence_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+make_flux_flux_equivalence_smoke)
  opm_add_test(make_flux_flux_equivalence_smoke
    DEPENDS
      flux_smoke_compare
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_make_flux_flux_equivalence_result_path}
      -f FLUX_DUMP_FLORES_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-make-flux-equivalence-smoke.sh ${PROJECT_BINARY_DIR}/bin/make_flux ${PROJECT_BINARY_DIR}/bin/flux_smoke_compare ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_FLORES_SMOKE.DATA FLUX_DUMP_FLORES_SMOKE flux 0"
  )

  set(_make_flux_flux_equivalence_mpi_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+make_flux_flux_equivalence_mpi_smoke)
  opm_add_test(make_flux_flux_equivalence_mpi_smoke
    DEPENDS
      flux_smoke_compare
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_make_flux_flux_equivalence_mpi_result_path}
      -f FLUX_DUMP_FLORES_SMOKE
      -n 2
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-make-flux-equivalence-smoke.sh ${PROJECT_BINARY_DIR}/bin/make_flux ${PROJECT_BINARY_DIR}/bin/flux_smoke_compare ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_FLORES_SMOKE.DATA FLUX_DUMP_FLORES_SMOKE flux 0"
  )

  set(_make_flux_pressure_sector_regression_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+make_flux_pressure_sector_regression_smoke)
  opm_add_test(make_flux_pressure_sector_regression_smoke
    DEPENDS
      flux_smoke_compare
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_make_flux_pressure_sector_regression_result_path}
      -f FLUX_DUMP_PSECTOR_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-make-flux-equivalence-smoke.sh ${PROJECT_BINARY_DIR}/bin/make_flux ${PROJECT_BINARY_DIR}/bin/flux_smoke_compare ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_PSECTOR_SMOKE.DATA FLUX_DUMP_PSECTOR_SMOKE pressure inline:box,8,8,3,22,13,3 --compare-arg=--ignore-times"
  )

  # The exterior relative permeabilities and capillary pressures make_flux
  # evaluates from the parent's saturation functions against the ones DUMPFLUX
  # reads off the converged state, on a deck whose rock changes at the sector
  # boundary.
  set(_make_flux_satfunc_equivalence_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+make_flux_satfunc_equivalence_smoke)
  opm_add_test(make_flux_satfunc_equivalence_smoke
    DEPENDS
      flux_smoke_compare
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_make_flux_satfunc_equivalence_result_path}
      -f FLUX_DUMP_SATFUNC_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-make-flux-equivalence-smoke.sh ${PROJECT_BINARY_DIR}/bin/make_flux ${PROJECT_BINARY_DIR}/bin/flux_smoke_compare ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_SATFUNC_SMOKE.DATA FLUX_DUMP_SATFUNC_SMOKE pressure inline:box,4,1,3,7,1,8 --compare-arg=--ignore-times"
    TEST_ARGS
      --initial-time-step-in-days=30
  )

  # ------------------------------------------------------------------------
  # USEFLUX: a sector run driven by the boundary data.

  set(_useflux_pressure_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_pressure_smoke)
  opm_add_test(useflux_pressure_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_pressure_result_path}
      -f FLUX_DUMP_PRESSURE_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-pressure-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_PRESSURE_SMOKE.DATA"
  )

  set(_useflux_flux_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_flux_smoke)
  opm_add_test(useflux_flux_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_flux_result_path}
      -f FLUX_DUMP_USEFLUX_FLUX_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-flux-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_FLUX_SMOKE.DATA"
  )

  set(_useflux_both_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_both_smoke)
  opm_add_test(useflux_both_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_both_result_path}
      -f FLUX_DUMP_USEFLUX_BOTH_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-both-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_BOTH_SMOKE.DATA"
  )

  set(_useflux_suffix_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_suffix_smoke)
  opm_add_test(useflux_suffix_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_suffix_result_path}
      -f FLUX_DUMP_USEFLUX_SUFFIX_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-suffix-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_SUFFIX_SMOKE.DATA"
  )

  set(_useflux_suffix2_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_suffix2_smoke)
  opm_add_test(useflux_suffix2_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_suffix2_result_path}
      -f FLUX_DUMP_USEFLUX_SUFFIX_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-suffix2-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_SUFFIX_SMOKE.DATA"
  )

  set(_useflux_missing_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_missing_smoke)
  opm_add_test(useflux_missing_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_missing_result_path}
      -f FLUX_DUMP_USEFLUX_SUFFIX_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-missing-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_MISSING_SMOKE.DATA"
  )

  set(_useflux_invalid_suffix_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_invalid_suffix_smoke)
  opm_add_test(useflux_invalid_suffix_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_invalid_suffix_result_path}
      -f FLUX_DUMP_USEFLUX_SUFFIX_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-invalid-suffix-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_INVALID_SUFFIX_SMOKE.DATA"
  )

  set(_useflux_missing_mpi_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_missing_mpi_smoke)
  opm_add_test(useflux_missing_mpi_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_missing_mpi_result_path}
      -f FLUX_DUMP_USEFLUX_SUFFIX_SMOKE
      -n 2
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-missing-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_MISSING_SMOKE.DATA 2"
  )

  set(_useflux_invalid_suffix_mpi_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_invalid_suffix_mpi_smoke)
  opm_add_test(useflux_invalid_suffix_mpi_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_invalid_suffix_mpi_result_path}
      -f FLUX_DUMP_USEFLUX_SUFFIX_SMOKE
      -n 2
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-invalid-suffix-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_INVALID_SUFFIX_SMOKE.DATA 2"
  )

  set(_useflux_mixed_suffix_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_mixed_suffix_smoke)
  opm_add_test(useflux_mixed_suffix_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_mixed_suffix_result_path}
      -f FLUX_DUMP_USEFLUX_SUFFIX_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-mixed-suffix-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_MIXED_SUFFIX_SMOKE.DATA"
  )

  set(_useflux_mixed_suffix_mpi_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_mixed_suffix_mpi_smoke)
  opm_add_test(useflux_mixed_suffix_mpi_smoke
    DEPENDS
      flux_smoke_shape_check
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_mixed_suffix_mpi_result_path}
      -f FLUX_DUMP_USEFLUX_SUFFIX_SMOKE
      -n 2
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-mixed-suffix-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_MIXED_SUFFIX_SMOKE.DATA 2"
  )

  set(_useflux_invalid_rate_nan_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_invalid_rate_nan_smoke)
  opm_add_test(useflux_invalid_rate_nan_smoke
    DEPENDS
      flux_smoke_shape_check
      flux_smoke_mutate
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_invalid_rate_nan_result_path}
      -f FLUX_DUMP_USEFLUX_FLUX_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-invalid-rate-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_BINARY_DIR}/bin/flux_smoke_mutate ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_FLUX_SMOKE.DATA nan"
  )

  set(_useflux_invalid_rate_inf_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_invalid_rate_inf_smoke)
  opm_add_test(useflux_invalid_rate_inf_smoke
    DEPENDS
      flux_smoke_shape_check
      flux_smoke_mutate
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_invalid_rate_inf_result_path}
      -f FLUX_DUMP_USEFLUX_FLUX_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-invalid-rate-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_BINARY_DIR}/bin/flux_smoke_mutate ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_FLUX_SMOKE.DATA inf"
  )

  set(_useflux_invalid_rate_nan_mpi_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_invalid_rate_nan_mpi_smoke)
  opm_add_test(useflux_invalid_rate_nan_mpi_smoke
    DEPENDS
      flux_smoke_shape_check
      flux_smoke_mutate
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_invalid_rate_nan_mpi_result_path}
      -f FLUX_DUMP_USEFLUX_FLUX_SMOKE
      -n 2
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-invalid-rate-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_BINARY_DIR}/bin/flux_smoke_mutate ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_FLUX_SMOKE.DATA nan 2"
  )

  set(_useflux_invalid_rate_inf_mpi_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+useflux_invalid_rate_inf_mpi_smoke)
  opm_add_test(useflux_invalid_rate_inf_mpi_smoke
    DEPENDS
      flux_smoke_shape_check
      flux_smoke_mutate
    EXE_TARGET
      flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_useflux_invalid_rate_inf_mpi_result_path}
      -f FLUX_DUMP_USEFLUX_FLUX_SMOKE
      -n 2
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-useflux-invalid-rate-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_BINARY_DIR}/bin/flux_smoke_mutate ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_FLUX_SMOKE.DATA inf 2"
  )

  if(OPM_SUMMARY_BIN)
    set(_flux_threshold_pressure_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_threshold_pressure_smoke)
    opm_add_test(flux_threshold_pressure_smoke
        EXE_TARGET flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/flux
        -r ${_flux_threshold_pressure_result_path}
        -f FLUX_DUMP_THPRES_SMOKE
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-flux-threshold-pressure.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_THPRES_SMOKE.DATA ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_THPRES_SMOKE.DATA ${PROJECT_BINARY_DIR}/bin/make_flux ${OPM_SUMMARY_BIN} ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_THPRES_CUTOUT_SMOKE.DATA"
    )

    set(_flux_fault_boundary_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_fault_boundary_smoke)
    opm_add_test(flux_fault_boundary_smoke
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/flux
        -r ${_flux_fault_boundary_result_path}
        -f FLUX_DUMP_FAULT_SMOKE
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-flux-fault-boundary-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_SUMMARY_BIN} ${PROJECT_BINARY_DIR}/bin/inspect_flux ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_FAULT_SMOKE.DATA ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_FAULT_SMOKE.DATA 1,1,5 2,1,5 I+"
    )

    set(_flux_fault_throw_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_fault_throw_smoke)
    opm_add_test(flux_fault_throw_smoke
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/flux
        -r ${_flux_fault_throw_result_path}
        -f FLUX_DUMP_THROW_SMOKE
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-flux-fault-throw.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_SUMMARY_BIN} ${PROJECT_BINARY_DIR}/bin/inspect_flux ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_THROW_SMOKE.DATA ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_THROW_SMOKE.DATA ${PROJECT_BINARY_DIR}/bin/make_flux"
    )

    set(_flux_boundary_interpolation_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_boundary_interpolation_smoke)
    opm_add_test(flux_boundary_interpolation_smoke
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_flux_boundary_interpolation_result_path}
        -f TWOWELL_PARENT_SPARSE_1D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-flux-boundary-interpolation.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_BINARY_DIR}/bin/make_flux ${OPM_SUMMARY_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_PARENT_SPARSE_1D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_USEFLUX_SPARSE_1D.DATA"
    )

    # A consumer must give the same answer however many ranks it runs on.
    # The NNC sector is used because it exercises both kinds of boundary
    # face, the axis-aligned ones and the one that arrives as a source term.
    set(_flux_parallel_consumer_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_parallel_consumer_smoke)
    opm_add_test(flux_parallel_consumer_smoke
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_flux_parallel_consumer_result_path}
        -f TWOWELL_DUMPFLUX_NNC_1D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-flux-parallel-consumer.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_SUMMARY_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_DUMPFLUX_NNC_1D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_USEFLUX_NNC_1D.DATA 2 WBHP:P1 FPR FOPT FGPT"
    )

    # Same seal, but on a K face, which a cell reaches through a different
    # slot in its boundary-face list because its I neighbours are inside the
    # sector. That is the arrangement the column case above cannot produce.
    set(_flux_kfault_boundary_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_kfault_boundary_smoke)
    opm_add_test(flux_kfault_boundary_smoke
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/flux
        -r ${_flux_kfault_boundary_result_path}
        -f FLUX_DUMP_KFAULT_SMOKE
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-flux-fault-boundary-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_SUMMARY_BIN} ${PROJECT_BINARY_DIR}/bin/inspect_flux ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_KFAULT_SMOKE.DATA ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_KFAULT_SMOKE.DATA 5,1,1 5,1,2 K+"
    )
  endif()

  if(OPM_COMPARE_ECL_BIN)
    set(_onewell_useflux_pressure_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+onewell_useflux_pressure_smoke)
    opm_add_test(onewell_useflux_pressure_smoke
      DEPENDS
        flux_smoke_shape_check
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_onewell_useflux_pressure_result_path}
        -f ONEWELL_DUMPFLUX_PRESSURE_1D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-onewell-useflux-pressure-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_COMPARE_ECL_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/ONEWELL_DUMPFLUX_PRESSURE_1D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/ONEWELL_USEFLUX_PRESSURE_1D.DATA"
    )

    set(_onewell_useflux_flux_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+onewell_useflux_flux_smoke)
    opm_add_test(onewell_useflux_flux_smoke
      DEPENDS
        flux_smoke_shape_check
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_onewell_useflux_flux_result_path}
        -f ONEWELL_DUMPFLUX_FLUX_1D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-onewell-useflux-flux-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_COMPARE_ECL_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/ONEWELL_DUMPFLUX_FLUX_1D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/ONEWELL_USEFLUX_FLUX_1D.DATA"
    )
  endif()

  # ------------------------------------------------------------------------
  # USEFLUX: the parent's wells, summary vectors and UDQs in a sector run.

  set(_flux_partial_regions_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_partial_regions_smoke)
  opm_add_test(flux_partial_regions_smoke
      EXE_TARGET flow_blackoil
    DRIVER_ARGS
      -i ${PROJECT_SOURCE_DIR}/tests/flux
      -r ${_flux_partial_regions_result_path}
      -f FLUX_DUMP_PARTIAL_REGIONS_SMOKE
      -p "bash ${PROJECT_SOURCE_DIR}/tests/check-flux-partial-regions.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_PARTIAL_REGIONS_SMOKE.DATA ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_PARTIAL_REGIONS_SMOKE.DATA 2"
  )

  if(OPM_SUMMARY_BIN)
    set(_make_flux_summary_cadence_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+make_flux_summary_cadence_smoke)
    opm_add_test(make_flux_summary_cadence_smoke
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_make_flux_summary_cadence_result_path}
        -f TWOWELL_PARENT_SPARSE_1D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-make-flux-summary-cadence.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_BINARY_DIR}/bin/make_flux ${PROJECT_BINARY_DIR}/bin/inspect_flux ${OPM_SUMMARY_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_PARENT_SPARSE_1D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_USEFLUX_SPARSE_1D.DATA"
    )

    set(_flux_cutout_wells_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_cutout_wells_smoke)
    opm_add_test(flux_cutout_wells_smoke
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_flux_cutout_wells_result_path}
        -f TWOWELL_CUTOUT_PARENT_1D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-flux-cutout-wells.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_BINARY_DIR}/bin/make_flux ${OPM_SUMMARY_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_CUTOUT_PARENT_1D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_CUTOUT_SECTOR_1D.DATA"
    )

    # Same again, for a sector whose group target is a UDA. The correction
    # that removes the absent wells' share of that target reads two things
    # that used to be answered rank by rank: which wells the sector holds,
    # and the uncorrected value, which only the I/O rank has. Both are now
    # settled across the ranks, which also means this path carries
    # collectives and has to be reached on every rank -- if it were not,
    # this test would hang rather than fail.
    set(_flux_uda_group_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+flux_uda_group_smoke)
    opm_add_test(flux_uda_group_smoke
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/flux
        -r ${_flux_uda_group_result_path}
        -f FLUX_DUMP_UDA_GROUP_SMOKE
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-flux-parallel-consumer.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_SUMMARY_BIN} ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_DUMP_UDA_GROUP_SMOKE.DATA ${PROJECT_SOURCE_DIR}/tests/flux/FLUX_USE_UDA_GROUP_SMOKE.DATA 2 GOPT:G1 WOPT:P3 FOPT"
    )
  endif()

  if(OPM_COMPARE_ECL_BIN)
    set(_twowell_useflux_flux_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+twowell_useflux_flux_smoke)
    opm_add_test(twowell_useflux_flux_smoke
      DEPENDS
        flux_smoke_shape_check
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_twowell_useflux_flux_result_path}
        -f TWOWELL_DUMPFLUX_FLUX_1D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-twowell-useflux-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_COMPARE_ECL_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_DUMPFLUX_FLUX_1D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_USEFLUX_FLUX_1D.DATA"
    )

    set(_twowell_useflux_grup_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+twowell_useflux_grup_smoke)
    opm_add_test(twowell_useflux_grup_smoke
      DEPENDS
        flux_smoke_shape_check
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_twowell_useflux_grup_result_path}
        -f TWOWELL_DUMPFLUX_GRUP_1D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-twowell-useflux-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_COMPARE_ECL_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_DUMPFLUX_GRUP_1D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_USEFLUX_GRUP_1D.DATA"
    )

    set(_twowell_useflux_live_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+twowell_useflux_live_smoke)
    opm_add_test(twowell_useflux_live_smoke
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_twowell_useflux_live_result_path}
        -f TWOWELL_DUMPFLUX_LIVE_1D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-twowell-useflux-live-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_COMPARE_ECL_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_DUMPFLUX_LIVE_1D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_USEFLUX_LIVE_1D.DATA"
    )

    set(_twowell_useflux_nnc_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+twowell_useflux_nnc_smoke)
    opm_add_test(twowell_useflux_nnc_smoke
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_twowell_useflux_nnc_result_path}
        -f TWOWELL_DUMPFLUX_NNC_1D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-twowell-useflux-nnc-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_COMPARE_ECL_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_DUMPFLUX_NNC_1D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_USEFLUX_NNC_1D.DATA"
    )

    set(_twowell_useflux_pressure_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+twowell_useflux_pressure_smoke)
    opm_add_test(twowell_useflux_pressure_smoke
      DEPENDS
        flux_smoke_shape_check
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_twowell_useflux_pressure_result_path}
        -f TWOWELL_DUMPFLUX_PRESSURE_1D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-twowell-useflux-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_COMPARE_ECL_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_DUMPFLUX_PRESSURE_1D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_USEFLUX_PRESSURE_1D.DATA"
    )

    set(_twowell_useflux_smry_fallback_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+twowell_useflux_smry_fallback_smoke)
    opm_add_test(twowell_useflux_smry_fallback_smoke
      DEPENDS
        flux_smoke_mutate
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_twowell_useflux_smry_fallback_result_path}
        -f TWOWELL_DUMPFLUX_FLUX_1D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-twowell-useflux-smry-fallback-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${PROJECT_BINARY_DIR}/bin/flux_smoke_mutate ${OPM_COMPARE_ECL_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_DUMPFLUX_FLUX_1D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_USEFLUX_FLUX_1D.DATA"
    )

    # The region map is dimensioned over every cell of the grid, so it also
    # assigns a region to inactive cells. This case has an inactive cell inside
    # the region, which is what separates the region from the set of cells that
    # actually take part in the flow.
    set(_twowell_useflux_inactive_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+twowell_useflux_inactive_smoke)
    opm_add_test(twowell_useflux_inactive_smoke
      DEPENDS
        flux_smoke_shape_check
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_twowell_useflux_inactive_result_path}
        -f TWOWELL_DUMPFLUX_INACTIVE_2D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-twowell-useflux-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_COMPARE_ECL_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_DUMPFLUX_INACTIVE_2D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_USEFLUX_INACTIVE_2D.DATA"
    )

    # Boundary data is written per time step by default; the report-step cadence
    # of earlier versions remains available behind a switch.
    set(_twowell_useflux_cadence_result_path ${PROJECT_BINARY_DIR}/tests/results/flux/flow_blackoil+twowell_useflux_cadence_smoke)
    opm_add_test(twowell_useflux_cadence_smoke
      DEPENDS
        flux_smoke_shape_check
      EXE_TARGET
        flow_blackoil
      DRIVER_ARGS
        -i ${PROJECT_SOURCE_DIR}/tests/data/sector
        -r ${_twowell_useflux_cadence_result_path}
        -f TWOWELL_DUMPFLUX_FLUX_1D
        -p "bash ${PROJECT_SOURCE_DIR}/tests/check-twowell-useflux-cadence-smoke.sh ${PROJECT_BINARY_DIR}/bin/flow_blackoil ${OPM_COMPARE_ECL_BIN} ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_DUMPFLUX_FLUX_1D.DATA ${PROJECT_SOURCE_DIR}/tests/data/sector/TWOWELL_USEFLUX_FLUX_1D.DATA"
    )
  endif()
endif()
