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

project(prisms_pf CXX)

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

file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" PRISMS_PF_VERSION LIMIT_COUNT 1)

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

configure_file(include/prismspf/config.h.in include/prismspf/config.h @ONLY)

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

# Start with finding MPI

# =========================================================
# Configure targets
# =========================================================

print_subsection("Configuring build targets")

# =========================================================
# Installation
# =========================================================

print_subsection("Configuring installation targets")

# =========================================================
# Testing
# =========================================================

print_subsection("Configuring PRISMS-PF tests")

# =========================================================
# Final steps to print some nice messages
# =========================================================
