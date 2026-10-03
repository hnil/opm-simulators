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
      flux_smoke_compare
    SOURCES
      tests/flux/flux_smoke_compare.cpp
    LIBRARIES
      opmcommon
  )

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
endif()
