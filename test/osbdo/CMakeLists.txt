# --------------------------------------------------------------------------- #
#    CMake file for the testers of the OSBDO examples                         #
#                                                                             #
#    The testers need a :MILPSolver both for the subproblems of the           #
#    BendersBFunction and for the reference LP or QP.                         #
# --------------------------------------------------------------------------- #

if (TARGET SMS++::MILPSolver)
    foreach (ex mcf sc ics fl)
        add_executable(osbdo_${ex}_test test_osbdo_${ex}.cpp)
        target_link_libraries(osbdo_${ex}_test PRIVATE
                              SMS++::BundleSolver
                              SMS++::MILPSolver)
    endforeach ()

    # the small instances only: those of the paper are for measurements
    add_test(NAME osbdo_mcf_test
             COMMAND osbdo_mcf_test mcf-20-100-3.txt BSPar.txt LPPar.txt
                                    LPPar_inner.txt
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    # the same with five re-optimizations after changes of the coupling
    add_test(NAME osbdo_mcf_seq_test
             COMMAND osbdo_mcf_test mcf-20-100-3.txt BSPar.txt LPPar.txt
                                    LPPar_inner.txt 0 1e-6 5 0.2 1
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME osbdo_sc_test
             COMMAND osbdo_sc_test sc-5.txt BSPar.txt LPPar.txt
                                   LPPar_inner.txt
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME osbdo_ics_test
             COMMAND osbdo_ics_test ics-20-30-4.txt BSPar.txt LPPar.txt
                                    LPPar_inner.txt
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME osbdo_fl_test
             COMMAND osbdo_fl_test fl-400-20-4.txt BSPar.txt
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif ()

# --------------------------------------------------------------------------- #

# ----- the label of this module's tests ------------------------------------ #

# The label lets the module CI select exactly the tests exercising this module
# with `ctest -L <module>` (see tests/cmake/TestLabels.cmake in the umbrella
# for the system-tests counterpart): without it these tests are built and
# never run.
get_property(_tests DIRECTORY PROPERTY TESTS)
if (_tests)
    set_tests_properties(${_tests} PROPERTIES LABELS "${modName}")
endif ()
