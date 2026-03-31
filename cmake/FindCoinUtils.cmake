# --------------------------------------------------------------------------- #
#    CMake find module for CoinUtils                                          #
#                                                                             #
#    This module finds CoinUtils include directories and libraries.           #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(CoinUtils [version] [EXACT] [REQUIRED])                 #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        CoinUtils_FOUND         - True if headers are found                  #
#        CoinUtils_INCLUDE_DIRS  - Include directories                        #
#        CoinUtils_LIBRARIES     - Libraries to be linked                     #
#        CoinUtils_VERSION       - Version number                             #
#                                                                             #
#    The search results are saved in these persistent cache entries:          #
#                                                                             #
#        CoinUtils_INCLUDE_DIR   - Directory containing headers               #
#        CoinUtils_LIBRARY       - The found library                          #
#        CoinUtils_DLL           - The found runtime DLL (Windows only)       #
#                                                                             #
#    This module can read a search path from the variable:                    #
#                                                                             #
#        CoinUtils_ROOT          - Preferred CoinUtils location               #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        Coin::CoinUtils                                                      #
#                                                                             #
#    This find module is provided because CoinUtils does not provide          #
#    a CMake configuration file on its own.                                   #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Requirements -------------------------------------------------------- #
find_package(BZip2 REQUIRED QUIET)

# Conda coin-or-utils package links MKL BLAS library
# https://github.com/conda-forge/coin-or-utils-feedstock/blob/main/recipe/build.sh#L12
if (WIN32 AND DEFINED ENV{CONDA_BUILD})
    find_library(MKL_RT_LIBRARY
            NAMES mkl_rt
            PATHS $ENV{LIBRARY_LIB}
            NO_DEFAULT_PATH
            DOC "MKL_RT library.")
endif ()

# Check if already in cache
if (WIN32)
    if (CoinUtils_INCLUDE_DIR AND CoinUtils_LIBRARY AND CoinUtils_DLL)
        set(CoinUtils_FOUND TRUE)
    endif ()
else ()
    if (CoinUtils_INCLUDE_DIR AND CoinUtils_LIBRARY)
        set(CoinUtils_FOUND TRUE)
    endif ()
endif ()

if (NOT CoinUtils_FOUND)

    # ----- Find the headers ------------------------------------------------ #
    find_path(CoinUtils_INCLUDE_DIR
            NAMES CoinUtilsConfig.h
            PATHS ${CoinUtils_ROOT}/include
            PATH_SUFFIXES coin coinutils/coin coin-or
            DOC "CoinUtils include directory.")

    # ----- Find the library ------------------------------------------------ #
    find_library(CoinUtils_LIBRARY
            NAMES CoinUtils
            PATHS ${CoinUtils_ROOT}/lib
            DOC "CoinUtils library.")

    # ----- Find the runtime DLL on Windows --------------------------------- #
    if (WIN32)
        find_file(CoinUtils_DLL
                NAMES CoinUtils.dll libCoinUtils.dll CoinUtils-0.dll
                PATHS ${CoinUtils_ROOT}/bin
                DOC "CoinUtils runtime DLL.")
    endif ()

    # ----- Parse the version ----------------------------------------------- #
    if (CoinUtils_INCLUDE_DIR)
        file(STRINGS
                "${CoinUtils_INCLUDE_DIR}/CoinUtilsConfig.h"
                _coinutils_version_lines REGEX "#define COINUTILS_VERSION_(MAJOR|MINOR|RELEASE)")

        string(REGEX REPLACE ".*COINUTILS_VERSION_MAJOR *\([0-9]*\).*" "\\1" _coinutils_version_major "${_coinutils_version_lines}")
        string(REGEX REPLACE ".*COINUTILS_VERSION_MINOR *\([0-9]*\).*" "\\1" _coinutils_version_minor "${_coinutils_version_lines}")
        string(REGEX REPLACE ".*COINUTILS_VERSION_RELEASE *\([0-9]*\).*" "\\1" _coinutils_version_release "${_coinutils_version_lines}")

        set(CoinUtils_VERSION "${_coinutils_version_major}.${_coinutils_version_minor}.${_coinutils_version_release}")
        unset(_coinutils_version_lines)
        unset(_coinutils_version_major)
        unset(_coinutils_version_minor)
        unset(_coinutils_version_release)
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    if (WIN32)
        find_package_handle_standard_args(
                CoinUtils
                REQUIRED_VARS CoinUtils_LIBRARY CoinUtils_DLL CoinUtils_INCLUDE_DIR
                VERSION_VAR CoinUtils_VERSION)
    else ()
        find_package_handle_standard_args(
                CoinUtils
                REQUIRED_VARS CoinUtils_LIBRARY CoinUtils_INCLUDE_DIR
                VERSION_VAR CoinUtils_VERSION)
    endif ()
endif ()

# ----- Export the target --------------------------------------------------- #
if (CoinUtils_FOUND)
    set(CoinUtils_INCLUDE_DIRS ${CoinUtils_INCLUDE_DIR})
    set(CoinUtils_LIBRARIES ${CoinUtils_LIBRARY})

    if (NOT TARGET Coin::CoinUtils)
        if (WIN32)
            add_library(Coin::CoinUtils SHARED IMPORTED)
            set_target_properties(
                    Coin::CoinUtils PROPERTIES
                    IMPORTED_IMPLIB "${CoinUtils_LIBRARY}"
                    IMPORTED_LOCATION "${CoinUtils_DLL}"
                    INTERFACE_INCLUDE_DIRECTORIES "${CoinUtils_INCLUDE_DIRS}")
        else ()
            add_library(Coin::CoinUtils UNKNOWN IMPORTED)
            set_target_properties(
                    Coin::CoinUtils PROPERTIES
                    IMPORTED_LOCATION "${CoinUtils_LIBRARY}"
                    INTERFACE_INCLUDE_DIRECTORIES "${CoinUtils_INCLUDE_DIRS}")
        endif ()

        target_link_libraries(Coin::CoinUtils INTERFACE "BZip2::BZip2")
        if (MKL_RT_LIBRARY)
            target_link_libraries(Coin::CoinUtils INTERFACE ${MKL_RT_LIBRARY})
        endif ()
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
if (WIN32)
    mark_as_advanced(CoinUtils_INCLUDE_DIR
            CoinUtils_LIBRARY
            CoinUtils_DLL
            CoinUtils_VERSION)
else ()
    mark_as_advanced(CoinUtils_INCLUDE_DIR
            CoinUtils_LIBRARY
            CoinUtils_VERSION)
endif ()

# --------------------------------------------------------------------------- #
