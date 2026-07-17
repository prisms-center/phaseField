##
#  CMake for the PRISMS-PF core library
##

# =========================================================
# CMake setup & loading functions
# =========================================================

# CMake config
message(STATUS "Using CMake ${CMAKE_VERSION}")
cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25)

# Prohibit in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(
    FATAL_ERROR
    "In-source build prohibited. \n"
    "Use a separate build directory with:\n  cmake -B build"
  )
endif()

# Set the build type to DebugRelease if unspecified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "DebugRelease" CACHE STRING "Build type" FORCE)
endif()

# Add DebugRelease, Debug, and Release to build types for
# interactive CMake
get_property(_types CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS)
list(
  APPEND _types
  DebugRelease
  Debug
  Release
)
list(REMOVE_DUPLICATES _types)
set_property(
  CACHE
    CMAKE_BUILD_TYPE
  PROPERTY
    STRINGS
      "${_types}"
)

# Project declaration, language, and version
project(prisms_pf LANGUAGES CXX VERSION 4.0.0)

# Export compile command for IDE lsp support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set the location of additional CMake modules
set(
  CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  ${PROJECT_SOURCE_DIR}/cmake
  ${PROJECT_SOURCE_DIR}/cmake/functions
)

include(prisms_pf_cmake_functions)

# =========================================================
# Configuration setup
# =========================================================

set(PRISMS_PF_VERSION ${PROJECT_VERSION})

print_subsection("Configuring PRISMS-PF v${PRISMS_PF_VERSION}")

include(prisms_pf_configuration)

message(STATUS "All good!")

# =========================================================
# Find dependencies
# =========================================================

print_subsection("Configuring external libraries")

# First, vendored dependencies
include(FetchContent)

FetchContent_Declare(
  libassert
  GIT_REPOSITORY https://github.com/jeremy-rifkin/libassert.git
  GIT_TAG v2.2.1
  SYSTEM
)

if(PRISMS_PF_UNIT_TESTS)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.13.0
    SYSTEM
  )
endif()

# Find MPI
include(deps/mpi)

# Find deal.II
include(deps/dealii)

# Find Caliper if we've enabled it
if(PRISMS_PF_WITH_CALIPER)
  include(deps/caliper)
endif()

# TODO: Add VTK

# FetchCotent MakeAvailable for the other vendored dependencies

# Okay, now we have to do a bunch of stuff because "DebugRelease" exists. It's
# nice, don't get me wrong, but it makes things way complicated, especially with
# vendored dependencies. Sure, we could just build the Release version of the
# dependencies but that defeats the point.
#
# To any future developers I hope you don't have to deal with this CMake :(
if(CMAKE_BUILD_TYPE STREQUAL "DebugRelease")
  # If we have DebugRelease as the build type, we want to build
  # libassert with Debug and Release targets to link to later
  include(ExternalProject)

  # Populate libassert so we can get the source directory
  FetchContent_Populate(libassert)

  # External projects for the debug and release libassert targets
  ExternalProject_Add(
    libassert_debug
    SOURCE_DIR ${libassert_SOURCE_DIR}
    PREFIX ${CMAKE_BINARY_DIR}/_deps/libassert-debug
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/libassert-debug/install
  )

  ExternalProject_Add(
    libassert_release
    SOURCE_DIR ${libassert_SOURCE_DIR}
    PREFIX ${CMAKE_BINARY_DIR}/_deps/libassert-release
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/libassert-release/install
  )

  # Create install directories
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/_deps/libassert-debug/install/include)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/_deps/libassert-release/install/include)

  # Add a custom library for libassert so we can link properly
  add_library(libassert_imported INTERFACE IMPORTED GLOBAL)
  add_dependencies(
    libassert_imported
    libassert_debug
    libassert_release
  )

  set_property(
    TARGET
      libassert_imported
    PROPERTY
      INTERFACE_LINK_LIBRARIES
        "$<$<CONFIG:Debug>:${CMAKE_BINARY_DIR}/_deps/libassert-debug/install/lib/libassert.a>"
        "$<$<CONFIG:Release>:${CMAKE_BINARY_DIR}/_deps/libassert-release/install/lib/libassert.a>"
  )
  set_property(
    TARGET
      libassert_imported
    PROPERTY
      INTERFACE_INCLUDE_DIRECTORIES
        "${CMAKE_BINARY_DIR}/_deps/libassert-debug/install/include"
        "${CMAKE_BINARY_DIR}/_deps/libassert-release/install/include"
  )

  add_library(libassert::assert ALIAS libassert_imported)

  # Custom install for libassert
  include(GNUInstallDirs)

  install(
    DIRECTORY
      ${CMAKE_BINARY_DIR}/_deps/libassert-debug/install/
    DESTINATION ${CMAKE_INSTALL_BINDIR}/../vendored/libassert-debug
  )

  install(
    DIRECTORY
      ${CMAKE_BINARY_DIR}/_deps/libassert-release/install/
    DESTINATION ${CMAKE_INSTALL_BINDIR}/../vendored/libassert
  )

  # Still make available other vendored packages (e.g. Catch2)
  # TODO: How should we have the build type here. Just whatever
  # the default is?
  list(REMOVE_ITEM PRISMS_PF_VENDORED_PACKAGES libassert)
  FetchContent_MakeAvailable(${PRISMS_PF_VENDORED_PACKAGES})
else()
  FetchContent_MakeAvailable(${PRISMS_PF_VENDORED_PACKAGES})
endif()

# =========================================================
# Configure targets
# =========================================================

print_subsection("Configuring build targets")
# TODO: Add status messages here

# Configure the config file
configure_file(include/prismspf/config.h.in include/prismspf/config.h @ONLY)

# Configure the template inst file and add the code generator script
add_subdirectory(cmake/scripts)
if(PRISMS_PF_ADDITIONAL_DEGREES)
  set(_additional_degrees "; 3; 4; 5; 6")
else()
  set(_additional_degrees "")
endif()
configure_file(cmake/templates.in cmake/templates @ONLY)

if(CMAKE_BUILD_TYPE STREQUAL "DebugRelease")
  add_library(prisms_pf_debug)
  add_library(prisms_pf)

  prisms_pf_configure_target(prisms_pf_debug "Debug")
  prisms_pf_configure_target(prisms_pf "Release")

  add_library(prisms_pf::prisms_pf_debug ALIAS prisms_pf_debug)
  add_library(prisms_pf::prisms_pf ALIAS prisms_pf)
else()
  add_library(prisms_pf)
  prisms_pf_configure_target(prisms_pf ${CMAKE_BUILD_TYPE})
  add_library(prisms_pf::prisms_pf ALIAS prisms_pf)
endif()

add_subdirectory(src)

# =========================================================
# Installation
# =========================================================

print_subsection("Configuring installation targets")
# TODO: Add status messages here

include(GNUInstallDirs)

install(
  TARGETS
    ${PRISMS_PF_TARGETS}
  EXPORT prisms_pf_targets
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
  FILE_SET HEADERS
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install the generated config.h separately since it's in the build tree
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/include/prismspf/config.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/prismspf
)

# Export the targets file so downstream can use find_package(PRISMS-PF)
install(
  EXPORT prisms_pf_targets
  FILE prisms_pf_targets.cmake
  NAMESPACE PRISMS-PF::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PRISMS-PF
)

# Generate and install the config file for find_package support
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/PRISMS-PFConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/PRISMS-PFConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PRISMS-PF
)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/PRISMS-PFConfigVersion.cmake
  VERSION ${PRISMS_PF_VERSION}
  COMPATIBILITY SameMinorVersion
)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/PRISMS-PFConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/PRISMS-PFConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PRISMS-PF
)
# Install executable utilities
install(
  PROGRAMS
    ${CMAKE_SOURCE_DIR}/contrib/postprocessing/visualize2D
  DESTINATION ${CMAKE_INSTALL_BINDIR}
)

export(
  TARGETS
    ${PRISMS_PF_TARGETS}
  NAMESPACE PRISMS-PF::
  FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/prisms_pf_targets.cmake"
)

# =========================================================
# Examples
# =========================================================

print_subsection("Configuring PRISMS-PF examples")
# TODO: Add status messages here

if(PRISMS_PF_EXAMPLES)
  add_subdirectory(applications)
endif()

# =========================================================
# Testing
# =========================================================

print_subsection("Configuring PRISMS-PF tests")
# TODO: Add status messages here

if(PRISMS_PF_UNIT_TESTS OR PRISMS_PF_REGRESSION_TESTS OR PRISMS_PF_PERFORMANCE_TESTS)
  enable_testing()
endif()

if(PRISMS_PF_UNIT_TESTS)
  add_subdirectory(tests/unit)
endif()
if(PRISMS_PF_REGRESSION_TESTS)
  #add_subdirectory(tests/regression)
endif()
if(PRISMS_PF_PERFORMANCE_TESTS)
  #add_subdirectory(tests/performance)
endif()

# =========================================================
# Documentation
# =========================================================

print_subsection("Configuring PRISMS-PF documentation")
# TODO: Add status message here

if(PRISMS_PF_DOCS)
  add_subdirectory(doc)
endif()

# =========================================================
# Final steps to print some nice messages
# =========================================================

file(READ ${CMAKE_BINARY_DIR}/summary.log PRISMS_PF_LOG_SUMMARY)
message("${PRISMS_PF_LOG_SUMMARY}")
