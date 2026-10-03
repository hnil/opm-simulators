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
endif()
