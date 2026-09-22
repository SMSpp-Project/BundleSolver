# --------------------------------------------------------------------------- #
#    CMake file for the tester of the large-scale functions of Karmitsa,      #
#    Bagirov and Makela                                                       #
#                                                                             #
#    The tester is linked with tnsunc.f, the file of test problems            #
#    distributed with the limited memory bundle method LMBM of Karmitsa,      #
#    which is not part of this repository: it is built only if KBM_TNSUNC    #
#    gives the path of that file.                                             #
# --------------------------------------------------------------------------- #

set(KBM_TNSUNC "" CACHE FILEPATH "path of tnsunc.f, the test problems of LMBM")

if (KBM_TNSUNC AND TARGET SMS++::MILPSolver)
    enable_language(Fortran)
    add_executable(kbm_test test_kbm.cpp ${KBM_TNSUNC})
    set_source_files_properties(${KBM_TNSUNC} PROPERTIES
                                COMPILE_OPTIONS "-std=legacy")
    target_link_libraries(kbm_test PRIVATE
                          SMS++::BundleSolver
                          SMS++::MILPSolver)
    if (BUILD_SHARED_LIBS AND UNIX AND (NOT APPLE))
        target_link_options(kbm_test PRIVATE "-Wl,--no-as-needed")
    endif ()
endif ()

# --------------------------------------------------------------------------- #
