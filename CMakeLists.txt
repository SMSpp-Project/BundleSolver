# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ BundleSolver                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.8...3.15)
project(BundleSolver VERSION 0.1.2 LANGUAGES CXX)

# ----- Settings ------------------------------------------------------------ #
include(CTest)

# ----- Requirements -------------------------------------------------------- #
find_package(SMSpp REQUIRED)
find_package(MILPSolver REQUIRED)
find_package(NDOFiOracle REQUIRED)

# ----- Library ------------------------------------------------------------- #
add_library(BundleSolver)
target_compile_features(BundleSolver PUBLIC cxx_std_17)

target_sources(BundleSolver PRIVATE BundleSolver.cpp)

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake), that must be different
# for local build and installation.
target_include_directories(
        BundleSolver PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/SMS++>)

target_link_libraries(BundleSolver PUBLIC
                      SMS++::SMSpp
                      SMS++::MILPSolver
                      NDOFiOracle::NDOFiOracle)

# Alias for testing
add_library(SMS++::BundleSolver ALIAS BundleSolver)

# ----- Install instructions ------------------------------------------------ #
# For installing the library and all its CMake modules on the system.
# Not required for local builds (see below).

include(GNUInstallDirs)

# Install library
install(TARGETS BundleSolver
        EXPORT BundleSolver-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install headers
install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/BundleSolver.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# Install target import file
install(EXPORT BundleSolver-targets
        FILE BundleSolverTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BundleSolver)

# Generate package version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/BundleSolverConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate package configuration file
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/BundleSolverConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/BundleSolverConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BundleSolver)

# Install package version and configuration files
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/BundleSolverConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/BundleSolverConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BundleSolver)

# ----- Export from the build tree ------------------------------------------ #
# This is needed so other projects can find the library in a local build.
# It does not affect installation.

if (CMAKE_EXPORT_PACKAGE_REGISTRY)

    export(EXPORT BundleSolver-targets
           FILE ${CMAKE_CURRENT_BINARY_DIR}/BundleSolverTargets.cmake
           NAMESPACE SMS++::)

    # Register package in user's package registry (~/.cmake/registry)
    export(PACKAGE BundleSolver)
endif ()

# ----- Subdirectories ------------------------------------------------------ #
if (BUILD_TESTING)
    add_subdirectory(Main)
endif()
