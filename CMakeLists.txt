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

include(ExternalProject)
include(GNUInstallDirs)

# Find libassert
include(deps/libassert)

# Find MPI
include(deps/mpi)

# Find deal.II
include(deps/dealii)

# Find VTK if we've enabled it
if(PRISMS_PF_WITH_VTK)
  include(deps/vtk)
endif()

# Find HDF5 if we've enabled it
if(PRISMS_PF_WITH_HDF5)
  include(deps/hdf5)
endif()

# Find Caliper if we've enabled it
if(PRISMS_PF_WITH_CALIPER)
  include(deps/caliper)
endif()

# Find Catch2 if we've enabled it
if(PRISMS_PF_UNIT_TESTS OR PRISMS_PF_REGRESSION_TESTS OR PRISMS_PF_PERFORMANCE_TESTS)
  include(deps/catch2)
endif()

# =========================================================
# Configure targets
# =========================================================

print_subsection("Configuring build targets")

# Configure the config file
configure_file(include/prismspf/config.h.in include/prismspf/config.h @ONLY)

# Configure the template inst file and add the code generator script
# TODO: Fix this to be not silly
add_subdirectory(cmake/scripts)
if(PRISMS_PF_ADDITIONAL_DEGREES)
  set(_additional_degrees "; 3; 4; 5; 6")
else()
  set(_additional_degrees "")
endif()
configure_file(cmake/templates.in cmake/templates @ONLY)

# TODO: Move this
# Create the compiler flags we care about
set(
  PRISMS_PF_CXX_FLAGS
  ${DEAL_II_CXX_FLAGS}
  ${CMAKE_CXX_FLAGS}
)
set(
  PRISMS_PF_CXX_FLAGS_DEBUG
  ${DEAL_II_CXX_FLAGS_DEBUG}
  ${CMAKE_CXX_FLAGS_DEBUG}
)
set(
  PRISMS_PF_CXX_FLAGS_RELEASE
  ${DEAL_II_CXX_FLAGS_RELEASE}
  ${CMAKE_CXX_FLAGS_RELEASE}
)

# Separate space delimited arguments into list
separate_arguments(PRISMS_PF_CXX_FLAGS_LIST UNIX_COMMAND "${PRISMS_PF_CXX_FLAGS}")
separate_arguments(
  PRISMS_PF_CXX_FLAGS_DEBUG_LIST
  UNIX_COMMAND
  "${PRISMS_PF_CXX_FLAGS_DEBUG}"
)
separate_arguments(
  PRISMS_PF_CXX_FLAGS_RELEASE_LIST
  UNIX_COMMAND
  "${PRISMS_PF_CXX_FLAGS_RELEASE}"
)

separate_arguments(DEAL_II_LINKER_FLAGS_LIST UNIX_COMMAND "${DEAL_II_LINKER_FLAGS}")
separate_arguments(
  DEAL_II_LINKER_FLAGS_DEBUG_LIST
  UNIX_COMMAND
  "${DEAL_II_LINKER_FLAGS_DEBUG}"
)
separate_arguments(
  DEAL_II_LINKER_FLAGS_RELEASE_LIST
  UNIX_COMMAND
  "${DEAL_II_LINKER_FLAGS_RELEASE}"
)

# TODO: Simplify this
if(CMAKE_BUILD_TYPE STREQUAL "DebugRelease")
  add_library(prisms_pf_debug)
  add_library(prisms_pf)
	add_dependencies(prisms_pf_debug prisms_pf)
  add_library(prisms_pf::prisms_pf_debug ALIAS prisms_pf_debug)
  add_library(prisms_pf::prisms_pf ALIAS prisms_pf)
else()
  add_library(prisms_pf)
  add_library(prisms_pf::prisms_pf ALIAS prisms_pf)
endif()

prisms_pf_configure_targets("${PRISMS_PF_TARGETS}")

add_subdirectory(src)

# =========================================================
# Installation
# =========================================================

print_subsection("Configuring installation targets")

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

# TODO: Install summary.log and detailed.log

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

# Export targets so you can link to build directory
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

if(PRISMS_PF_EXAMPLES)
  add_subdirectory(applications)
endif()

# =========================================================
# Testing
# =========================================================

print_subsection("Configuring PRISMS-PF tests")

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

if(PRISMS_PF_DOCS)
  add_subdirectory(doc)
endif()

# =========================================================
# Final steps to print some nice messages
# =========================================================

# TODO: Handle this
# file(READ ${CMAKE_BINARY_DIR}/summary.log PRISMS_PF_LOG_SUMMARY)
# message("${PRISMS_PF_LOG_SUMMARY}")
