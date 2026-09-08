# --------------------------------------------------------------------------- #
#    CMake file for BundleSolver/test                                         #
#                                                                             #
#    This is the internal test of the BundleSolver module: a small,           #
#    self-contained test of BundleSolverML that depends only on the SMS++     #
#    core and BundleSolver itself. Since BundleSolverML is only built into    #
#    BundleSolver when Torch is available, the test is added only if Torch    #
#    is found.                                                                #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

find_package(Torch QUIET HINTS ${Torch_ROOT})

if (Torch_FOUND)

    # Since we are using the block/solver factory, objects from linked
    # libraries may wrongly appear unused, and by default the linker does not
    # include them, so we have to force the linking.
    if (BUILD_SHARED_LIBS)
        if (UNIX AND (NOT APPLE))
            add_link_options("-Wl,--no-as-needed")
        endif ()
    else ()
        if (MSVC)
            add_link_options("/WHOLEARCHIVE" "/FORCE:MULTIPLE")
        else () # Unix
            if (APPLE)
                add_link_options("-Wl,-all_load")
            else ()
                add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
            endif ()
        endif ()
    endif ()

    # ----- BundleSolverML_test --------------------------------------------- #

    add_executable(${modName}ML_test test.cpp)
    target_link_libraries(${modName}ML_test PRIVATE
                          SMS++::BundleSolver
                          ${TORCH_LIBRARIES})

    add_executable(LSTMBlock_test LSTMBlock_test.cpp)
    target_link_libraries(LSTMBlock_test PRIVATE
                          SMS++::BundleSolver
                          ${TORCH_LIBRARIES})

    add_executable(BundleSolverMLCheckpoint_test BundleSolverMLCheckpoint_test.cpp)
    target_link_libraries(BundleSolverMLCheckpoint_test PRIVATE
                          SMS++::BundleSolver
                          ${TORCH_LIBRARIES})

    if (WIN32 AND BUILD_SHARED_LIBS AND NOT DEFINED ENV{CONDA_BUILD})
        # On Windows with shared libraries enabled copies all runtime DLL
        # dependencies next to the executable in the build directory
        # (POST_BUILD), using CMake's TARGET_RUNTIME_DLLS generator expression
        add_custom_command(TARGET ${modName}ML_test POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_RUNTIME_DLLS:${modName}ML_test>
                $<TARGET_FILE_DIR:${modName}ML_test>
                COMMAND_EXPAND_LISTS
                COMMENT "Copying runtime DLLs next to ${modName}ML_test.exe")
    endif ()

    # the test reads BSPar.txt and BSPar-ML.txt from the working directory
    add_test(NAME ${modName}ML_test
             COMMAND ${modName}ML_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

endif ()

# --------------------------------------------------------------------------- #

# ----- the quadratic 0-th component ---------------------------------------- #
# Neither of these needs Torch: the first checks the arithmetic the master
# absorbs the component with, the second solves the same regularised risk
# twice, with the bundle and as one quadratic program, and compares. The
# linker has to be told to keep the libraries whose objects only register
# themselves into the factories
add_executable(zeroth_quadratic_test test_zeroth_quadratic.cpp)
target_link_libraries(zeroth_quadratic_test PRIVATE SMS++::BundleSolver)
add_test(NAME zeroth_quadratic_test COMMAND zeroth_quadratic_test)

if (TARGET SMS++::MILPSolver)
    add_executable(zeroth_quadratic_e2e_test test_zeroth_quadratic_e2e.cpp)
    target_link_libraries(zeroth_quadratic_e2e_test PRIVATE
                          SMS++::BundleSolver
                          SMS++::MILPSolver)
    if (BUILD_SHARED_LIBS AND UNIX AND (NOT APPLE))
        target_link_options(zeroth_quadratic_e2e_test PRIVATE
                            "-Wl,--no-as-needed")
    endif ()
    # the test reads BSPar-rho.txt and MPBCfg-rho.txt from the working directory
    add_test(NAME zeroth_quadratic_e2e_test
             COMMAND zeroth_quadratic_e2e_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif ()
