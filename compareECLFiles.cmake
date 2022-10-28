# This script manages the addition of tests.
# The tests are orchestrated by a shell script,
# configured using opm_set_test_driver()
# and then the appropriate helper macro is called to
# register the ctest entry through the opm_add_test macro.
# Information such as the binary to call and test tolerances
# are passed from the build system to the driver script through
# command line parameters. See the opm_add_test() documentation for
# details on the parameters passed to the macro.

# Define some paths
set(BASE_RESULT_PATH ${PROJECT_BINARY_DIR}/tests/results)

###########################################################################
# TEST: runSim
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class simply runs a simulation.
function(add_test_runSimulator)
  set(oneValueArgs CASENAME FILENAME SIMULATOR DIR DIR_PREFIX PROCS)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  set(RESULT_PATH ${BASE_RESULT_PATH}${PARAM_DIR_PREFIX}/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${PARAM_TEST_ARGS})
  opm_add_test(runSimulator/${PARAM_CASENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                           -r ${RESULT_PATH}
                           -b ${PROJECT_BINARY_DIR}/bin
                           -f ${PARAM_FILENAME}
                           -n ${PARAM_PROCS}
               TEST_ARGS ${TEST_ARGS}
               CONFIGURATION extra)
endfunction()

###########################################################################
# TEST: compareECLFiles
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares output from a simulation to reference files.
function(add_test_compareECLFiles)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR DIR_PREFIX PREFIX RESTART_STEP RESTART_SCHED)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  if(NOT PARAM_PREFIX)
    set(PARAM_PREFIX compareECLFiles)
  endif()

  set(RESULT_PATH ${BASE_RESULT_PATH}${PARAM_DIR_PREFIX}/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${PARAM_TEST_ARGS})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c ${COMPARE_ECL_COMMAND}
                  -d ${RST_DECK_COMMAND})
   if(PARAM_RESTART_STEP)
     list(APPEND DRIVER_ARGS -s ${PARAM_RESTART_STEP})
   endif()
  if(PARAM_RESTART_SCHED)
   list(APPEND DRIVER_ARGS -h ${PARAM_RESTART_SCHED})
  endif()
  opm_add_test(${PARAM_PREFIX}_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${TEST_ARGS})
endfunction()

###########################################################################
# TEST: add_test_compare_restarted_simulation
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares the output from a restarted simulation
#     to that of a non-restarted simulation.
function(add_test_compare_restarted_simulation)
  set(oneValueArgs CASENAME FILENAME SIMULATOR TEST_NAME ABS_TOL REL_TOL DIR RESTART_STEP)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  if (PARAM_TEST_NAME)
    set(TEST_NAME ${PARAM_TEST_NAME})
  else()
    set(TEST_NAME compareRestartedSim_${PARAM_SIMULATOR}+${PARAM_FILENAME})
  endif()

  set(RESULT_PATH ${BASE_RESULT_PATH}/restart/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  opm_add_test(${TEST_NAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                           -r ${RESULT_PATH}
                           -b ${PROJECT_BINARY_DIR}/bin
                           -f ${PARAM_FILENAME}
                           -a ${PARAM_ABS_TOL}
                           -t ${PARAM_REL_TOL}
                           -c ${COMPARE_ECL_COMMAND}
                           -d ${RST_DECK_COMMAND}
                           -s ${PARAM_RESTART_STEP}
               TEST_ARGS ${PARAM_TEST_ARGS})
endfunction()

###########################################################################
# TEST: add_test_compare_parallel_simulation
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares the output from a parallel simulation
#     to the output from the serial instance of the same model.
function(add_test_compare_parallel_simulation)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()

  set(RESULT_PATH ${BASE_RESULT_PATH}/parallel/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${OPM_TESTS_ROOT}/${PARAM_DIR}/${PARAM_FILENAME} ${PARAM_TEST_ARGS})

  # Add test that runs flow_mpi and outputs the results to file
  opm_add_test(compareParallelSim_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                           -r ${RESULT_PATH}
                           -b ${PROJECT_BINARY_DIR}/bin
                           -f ${PARAM_FILENAME}
                           -a ${PARAM_ABS_TOL}
                           -t ${PARAM_REL_TOL}
                           -c ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${TEST_ARGS})
  set_tests_properties(compareParallelSim_${PARAM_SIMULATOR}+${PARAM_FILENAME}
                       PROPERTIES RUN_SERIAL 1)
endfunction()


###########################################################################
# TEST: add_test_compare_parallel_restarted_simulation
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares the output from a restarted parallel simulation
#     to that of a non-restarted parallel simulation.
function(add_test_compare_parallel_restarted_simulation)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR MPI_PROCS RESTART_STEP)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  set(RESULT_PATH ${BASE_RESULT_PATH}/parallelRestart/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c ${COMPARE_ECL_COMMAND}
                  -s ${PARAM_RESTART_STEP}
                  -d ${RST_DECK_COMMAND})
  if(PARAM_MPI_PROCS)
    list(APPEND DRIVER_ARGS -n ${PARAM_MPI_PROCS})
  endif()

  opm_add_test(compareParallelRestartedSim_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${PARAM_TEST_ARGS})
  set_tests_properties(compareParallelRestartedSim_${PARAM_SIMULATOR}+${PARAM_FILENAME}
                       PROPERTIES RUN_SERIAL 1)
endfunction()


###########################################################################
# TEST: add_test_split_comm
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares the output from a parallel simulation
#     to that of a parallel simulation running with a custom communicator.
function(add_test_split_comm)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  set(RESULT_PATH ${BASE_RESULT_PATH}/parallelSplitComm/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c ${COMPARE_ECL_COMMAND})
  if(PARAM_MPI_PROCS)
    list(APPEND DRIVER_ARGS -n ${PARAM_MPI_PROCS})
  endif()

  opm_add_test(compareParallelSplitComm_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${PARAM_TEST_ARGS})
  set_tests_properties(compareParallelSplitComm_${PARAM_SIMULATOR}+${PARAM_FILENAME}
                       PROPERTIES RUN_SERIAL 1)
endfunction()


###########################################################################


if(NOT TARGET test-suite)
  add_custom_target(test-suite)
endif()

# Simple execution tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-test.sh "")
add_test_runSimulator(CASENAME norne
                      FILENAME NORNE_ATW2013
                      SIMULATOR flow
                      PROCS 1)

add_test_runSimulator(CASENAME norne_parallel
                      FILENAME NORNE_ATW2013
                      SIMULATOR flow
                      DIR norne
                      PROCS 4)

include (${CMAKE_CURRENT_SOURCE_DIR}/regressionTests.cmake)
include (${CMAKE_CURRENT_SOURCE_DIR}/restartTests.cmake)

#set(cpr_abs_tol ${abs_tol})
set(cpr_abs_tol 1e99)
set(cpr_rel_tol ${rel_tol})                       
set(LINSOLVERS cpr_quasiimpes cpr_trueimpes cpr_trueimpesanalytic cpr cprw)
foreach(LINSOLVER IN LISTS LINSOLVERS)
   add_test_compareECLFiles(CASENAME spe12
                          FILENAME SPE1CASE2
                          SIMULATOR flow
                          PREFIX ${LINSOLVER}
                          ABS_TOL ${cpr_abs_tol}
                          REL_TOL ${cpr_rel_tol}
                          DIR spe1
                          TEST_ARGS --linsolver=${LINSOLVER})

endforeach()
add_test_compareECLFiles(CASENAME spe12
                         FILENAME SPE1CASE2
                         SIMULATOR flow
                         PREFIX cpr_addwell
                         ABS_TOL ${cpr_abs_tol}
                         REL_TOL ${cpr_rel_tol}
                         DIR spe1
                         TEST_ARGS --linsolver=cpr --matrix-add-well-contributions=true)
                       
# Restart tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-restart-regressionTest.sh "")
# Cruder tolerances for the restarted tests
set(abs_tol_restart 2e-1)
set(rel_tol_restart 4e-4)

add_test_compare_restarted_simulation(CASENAME spe1
                                      FILENAME SPE1CASE2_ACTNUM
                                      SIMULATOR flow
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      RESTART_STEP 6
                                      TEST_ARGS --sched-restart=false)

add_test_compare_restarted_simulation(CASENAME spe9
                                      FILENAME SPE9_CP_SHORT
                                      SIMULATOR flow
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      RESTART_STEP 15
                                      TEST_ARGS --sched-restart=false  --tolerance-mb=1e-7)

add_test_compare_restarted_simulation(CASENAME ctaquifer_2d_oilwater
                                      FILENAME 2D_OW_CTAQUIFER
                                      SIMULATOR flow
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      DIR aquifer-oilwater
                                      RESTART_STEP 15
                                      TEST_ARGS --sched-restart=true)

add_test_compare_restarted_simulation(CASENAME fetkovich_2d
                                      FILENAME 2D_FETKOVICHAQUIFER
                                      SIMULATOR flow
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      RESTART_STEP 30
                                      DIR aquifer-fetkovich
                                      TEST_ARGS --sched-restart=true)

add_test_compare_restarted_simulation(CASENAME numerical_aquifer_3d_1aqu
                                      FILENAME 3D_1AQU_3CELLS
                                      SIMULATOR flow
                                      ABS_TOL 0.4
                                      REL_TOL 4.0e-3
                                      RESTART_STEP 3
                                      DIR aquifer-num
                                      TEST_ARGS --sched-restart=true --enable-tuning=true)

add_test_compare_restarted_simulation(CASENAME numerical_aquifer_3d_2aqu
                                      FILENAME 3D_2AQU_NUM
                                      SIMULATOR flow
                                      ABS_TOL 0.4
                                      REL_TOL 4.0e-3
                                      RESTART_STEP 3
                                      DIR aquifer-num
                                      TEST_ARGS --sched-restart=true --enable-tuning=true)

# The dynamic MSW data is not written to /read from the restart file
# We therefore accept significant deviation in the results.
# Note also that we use --sched-restart=true since some necessary
# MSW info is still lacking in the restart file.
set(abs_tol_restart_msw 2e2)
set(rel_tol_restart_msw 1e-3)

add_test_compare_restarted_simulation(CASENAME msw_3d_hfa
                                      FILENAME 3D_MSW
                                      SIMULATOR flow
                                      ABS_TOL ${abs_tol_restart_msw}
                                      REL_TOL ${rel_tol_restart_msw}
                                      RESTART_STEP 10
                                      TEST_ARGS --enable-adaptive-time-stepping=false --sched-restart=true)


# Basic restart tests which only compare the summary output, this test driver should
# only be used in situations where it is challenging to get agreement in the restart file.
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-summary-restart-regressionTest.sh "")

add_test_compare_restarted_simulation(CASENAME spe1
                                      FILENAME SPE1CASE2_ACTNUM
                                      SIMULATOR flow
                                      TEST_NAME restart_spe1_summary
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      RESTART_STEP 6
                                      TEST_ARGS --sched-restart=false)

>>>>>>> origin/cpr_analytic

# PORV test
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-porv-acceptanceTest.sh "")
add_test_compareECLFiles(CASENAME norne
                         FILENAME NORNE_ATW2013
                         SIMULATOR flow
                         ABS_TOL 1e-5
                         REL_TOL 1e-8
                         PREFIX comparePORV
                         DIR_PREFIX /porv)

# Init tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-init-regressionTest.sh "")

add_test_compareECLFiles(CASENAME norne
                         FILENAME NORNE_ATW2013
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         PREFIX compareECLInitFiles
                         DIR_PREFIX /init)

# This is not a proper regression test; the test will load a norne case prepared
# for restart and run one single timestep - of length one day. The results are not
# verified in any way.
add_test(NAME NORNE_RESTART
         COMMAND flow --output-dir=${BASE_RESULT_PATH}/norne-restart ${OPM_TESTS_ROOT}/norne/NORNE_ATW2013_RESTART.DATA)


if(MPI_FOUND)
  include (${CMAKE_CURRENT_SOURCE_DIR}/parallelRestartTests.cmake)

  # Single test to verify that we treat custom communicators correctly.
  opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-split-comm-test.sh "")
  add_test_split_comm(CASENAME spe1
                      FILENAME SPE1CASE2
                      SIMULATOR flow
                      ABS_TOL 0.0
                      REL_TOL 0.0)

  include (${CMAKE_CURRENT_SOURCE_DIR}/parallelTests.cmake)
endif()


if(OPM_TESTS_ROOT)
  add_custom_target(update_data
                    COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target check --parallel || exit 0
                    COMMAND ${CMAKE_COMMAND} -E env "REASON=Local\ data\ update" ${PROJECT_SOURCE_DIR}/tests/update_reference_data.sh ${OPM_TESTS_ROOT} ${PROJECT_BINARY_DIR}
                    COMMENT "Updating reference data")
endif()
