##
#  CMake for the PRISMS-PF core library
##

# =========================================================
# CMake setup & loading functions
# =========================================================

# CMake config
message(STATUS "Using CMake ${CMAKE_VERSION}")
cmake_minimum_required(VERSION 3.23)
cmake_policy(VERSION 3.23)

# Prohibit in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source build prohibited.")
endif()

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
version_decomposition(${PRISMS_PF_VERSION})

print_subsection("Configuring PRISMS-PF v${PRISMS_PF_VERSION}")

prisms_pf_git_version()

print_variables(
  ""
  PRISMS_PF_GIT_BRANCH
  PRISMS_PF_GIT_TAG
  PRISMS_PF_GIT_REVISION
  PRISMS_PF_GIT_SHORTREV
  PRISMS_PF_GIT_TIMESTAMP
  ""
)

# TODO: Refactor this somewhere else
# Configuration options
option(PRISMS_PF_AUTODETECTION "Autodetection of PRISMS-PF dependencies." ON)

option(PRISMS_PF_UNIT_TESTS "Whether to build the unit tests or not." OFF)
option(PRISMS_PF_REGRESSION_TESTS "Whether to build the regression tests or not." OFF)
option(PRISMS_PF_PERFORMANCE_TESTS "Whether to build the performance tests or not." OFF)
option(PRISMS_PF_EXAMPLES "Whether to build examples or not." OFF)

# NOTE: This relies on the fact that deal.II is also compiled with 64-bit indices
option(PRISMS_PF_64BIT_INDICES "Whether to use 64-bit indices for large simulations" OFF)
# NOTE: This can only go higher than 6 if deal.II is compiled with >6 for the element degree.
option(PRISMS_PF_ADDITIONAL_DEGREES "Whether to compile with element degrees of 3+." OFF)
# NOTE: When on deal.II will try to use the max number of threads with respect to the number
# of MPI processes.
option(PRISMS_PF_THREADS "Whether to use threading." ON)

option(
  PRISMS_PF_WITH_VTK
  "Whether the user wants to compile PRISMS-PF with vtk, or not"
  OFF
)
option(
  PRISMS_PF_WITH_CALIPER
  "Whether the user wants to compile PRISMS-PF with the profiling code Caliper, or not."
  OFF
)

print_variables(
  ""
  PRISMS_PF_AUTODETECTION
  ""
  PRISMS_PF_UNIT_TESTS
  PRISMS_PF_REGRESSION_TESTS
  PRISMS_PF_PERFORMANCE_TESTS
  PRISMS_PF_EXAMPLES
  ""
  PRISMS_PF_64BIT_INDICES
  PRISMS_PF_ADDITIONAL_DEGREES
  PRISMS_PF_THREADS
  ""
  PRISMS_PF_WITH_VTK
  PRISMS_PF_WITH_CALIPER
  ""
)

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
)

if(PRISMS_PF_UNIT_TESTS)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.13.0
  )
endif()

# TODO: Refactor this
# Find MPI
# TODO: We should see if there's a minimum version we need to specify
find_package(MPI REQUIRED COMPONENTS CXX)

if(NOT MPI_CXX_FOUND)
  message(
    FATAL_ERROR
    "MPI CXX component not found. Please make sure your MPI installation includes C++ support."
  )
endif()
if(NOT MPI_CXX_COMPILER)
  message(
    FATAL_ERROR
    "MPI CXX compiler wrapper not found. Please make sure you have `mpicxx`"
  )
endif()

# Find the MPI CXX version if not already populated
if(NOT MPI_CXX_VERSION)
  file(READ "${MPI_CXX_INCLUDE_DIRS}/mpi.h" _mpi_header_content)

  string(REGEX MATCH "#define MPI_MAJOR_VERSION ([0-9]+)" _ "${_mpi_header_content}")
  set(_mpi_major ${CMAKE_MATCH_1})

  string(REGEX MATCH "#define MPI_MINOR_VERSION ([0-9]+)" _ "${_mpi_header_content}")
  set(_mpi_major ${CMAKE_MATCH_1})

  set(MPI_CXX_VERSION "${_mpi_major}.${_mpi_minor}")
endif()

# TODO: Refactor this

# Once we've found MPI, we need to find deal.II. Since deal.II likes to do things
# a little differently, provide some directory hints.
find_package(
  deal.II
  9.6.0
  REQUIRED
  HINTS
    ${DEAL_II_DIR}
    $ENV{DEAL_II_DIR}
)

# Again we run into some more of the uniqueness of deal.II with the way they set
# up build types. Their default build type is `DebugRelease`, which builds both
# Debug and Release targets. As such, our project much choose what deal.II target
# to link to when the user compiles with Debug and Release.
#
# When deal.II is built with only Debug or only Release, we link against for all
# PRISMS-PF builds. For example, if deal.II is built in Release, PRISMS-PF in
# Debug will only contain assertions and debugging for our library. In most cases,
# this is sufficient to catch most issues.
#
# When deal.II is build with DebugRelease, we link to the target that corresponds
# to the PRISMS-PF build type. This is recommended for developers.

# print_variable(DEAL_II_BUILD_TYPE)

set(PRISMS_PF_VENDORED_PACKAGES libassert)
if(PRISMS_PF_UNIT_TESTS OR PRISMS_PF_REGRESSION_TESTS OR PRISMS_PF_PERFORMANCE_TESTS)
  list(APPEND PRISMS_PF_VENDORED_PACKAGES Catch2)
endif()

FetchContent_MakeAvailable(${PRISMS_PF_VENDORED_PACKAGES})

# =========================================================
# Configure targets
# =========================================================

print_subsection("Configuring build targets")

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

# Declare the library
add_library(prisms_pf)

# Set the minimum C++ standard. This must be public so
# our applications also obey this.
target_compile_features(prisms_pf PUBLIC cxx_std_20)
set_target_properties(
  prisms_pf
  PROPERTIES
    CXX_EXTENSIONS
      OFF
)

# Add each modules in the src folder
add_subdirectory(src)

# Include the config to the library
target_include_directories(
  prisms_pf
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Link dependencies
target_link_libraries(
  prisms_pf
  PUBLIC
    MPI::MPI_CXX
    libassert::assert
    dealii::dealii
)

# =========================================================
# Installation
# =========================================================

print_subsection("Configuring installation targets")

include(GNUInstallDirs)

install(
  TARGETS
    prisms_pf
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

# Export the targets file so downstream can use find_package(prisms_pf)
install(
  EXPORT prisms_pf_targets
  FILE prisms_pf_targets.cmake
  NAMESPACE prisms_pf::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prisms_pf
)

# Generate and install the config file for find_package support
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/prisms_pf_config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/prisms_pf_config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prisms_pf
)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/prisms_pf_config_version.cmake
  VERSION ${PRISMS_PF_VERSION}
  COMPATIBILITY SameMinorVersion
)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/prisms_pf_config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/prisms_pf_config_version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prisms_pf
)

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
  add_subdirectory(tests/regression)
endif()
if(PRISMS_PF_PERFORMANCE_TESTS)
  add_subdirectory(tests/performance)
endif()

# =========================================================
# Final steps to print some nice messages
# =========================================================
