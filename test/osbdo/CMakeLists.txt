# --------------------------------------------------------------------------- #
#    CMake file for the multicommodity tester of the OSBDO instances          #
#                                                                             #
#    The tester needs a :MILPSolver both for the subproblems of the           #
#    BendersBFunction and for the reference LP.                               #
# --------------------------------------------------------------------------- #

if (TARGET SMS++::MILPSolver)
    add_executable(osbdo_mcf_test test_osbdo_mcf.cpp)
    target_link_libraries(osbdo_mcf_test PRIVATE
                          SMS++::BundleSolver
                          SMS++::MILPSolver)
    if (BUILD_SHARED_LIBS AND UNIX AND (NOT APPLE))
        target_link_options(osbdo_mcf_test PRIVATE "-Wl,--no-as-needed")
    endif ()
    # the small instance only: the one of the paper is for measurements
    add_test(NAME osbdo_mcf_test
             COMMAND osbdo_mcf_test mcf-20-100-3.txt BSPar.txt LPPar.txt
                                    LPPar_inner.txt
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif ()

# --------------------------------------------------------------------------- #
