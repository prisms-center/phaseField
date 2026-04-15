##
#  CMake for the PRISMS-PF core library
##

# =========================================================
# CMake setup & loading functions
# =========================================================

# CMake config
message(STATUS "Using CMake ${CMAKE_VERSION}")
cmake_minimum_required(VERSION 3.28)
cmake_policy(VERSION 3.28)

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
option(PRISMS_PF_DOCS "Whether to build documentation or not." OFF)
option(PRISMS_PF_CLANG_TIDY "Whether to run clang-tidy during the build stage." OFF)

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

if(PRISMS_PF_CLANG_TIDY)
  # NOTE: I've tested clang-tidy-18 and found that it fails. clang-tidy-20 works.
  find_program(
    CLANG_TIDY_TOOL
    NAMES
      clang-tidy
      clang-tidy-20
      clang-tidy-19
  )
  execute_process(
    COMMAND
      ${CLANG_TIDY_TOOL} --version
    OUTPUT_VARIABLE CLANG_TIDY_TOOL_VERSION
    ERROR_VARIABLE CLANG_TIDY_TOOL_VERSION
  )
  string(
    REGEX MATCH
    "[0-9]+(\\.[0-9]+)+"
    CLANG_TIDY_TOOL_VERSION
    "${CLANG_TIDY_TOOL_VERSION}"
  )
  message(STATUS "clang-tidy ${CLANG_TIDY_TOOL_VERSION} (${CLANG_TIDY_TOOL})")
  # NOTE: To avoid linting on the vendored source files we set CXX_CLANG_TIDY
  # as a property of our targets
endif()

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

# Print some information about deal.II configuration
print_variables(
  ""
  DEAL_II_VERSION
  DEAL_II_BUILD_TYPE
  DEAL_II_CXX_COMPILER
  DEAL_II_C_COMPILER
  ""
  DEAL_II_WITH_THREADS
  DEAL_II_WITH_64BIT_INDICES
  DEAL_II_WITH_BOOST
  DEAL_II_BOOST_VERSION
  DEAL_II_WITH_KOKKOS
  DEAL_II_KOKKOS_VERSION
  DEAL_II_FEATURE_KOKKOS_BUNDLED_CONFIGURED
  DEAL_II_WITH_LAPACK
  DEAL_II_LAPACK_WITH_MKL
  DEAL_II_WITH_MPI
  DEAL_II_MPI_VERSION
  DEAL_II_WITH_P4EST
  DEAL_II_P4EST_VERSION
  DEAL_II_P4EST_WITH_VTK_BINARY
  DEAL_II_WITH_TASKFLOW
  DEAL_II_TASKFLOW_VERSION
  DEAL_II_FEATURE_TASKFLOW_BUNDLED_CONFIGURED
  DEAL_II_WITH_TBB
  DEAL_II_WITH_VTK
  DEAL_II_WITH_ZLIB
  DEAL_II_ZLIB_VERSION
  ""
  DEAL_II_CXX_FLAGS
  DEAL_II_CXX_FLAGS_DEBUG
  DEAL_II_CXX_FLAGS_RELEASE
  ""
  DEAL_II_WARNING_FLAGS
  ""
  DEAL_II_LINKER_FLAGS
  DEAL_II_LINKER_FLAGS_DEBUG
  DEAL_II_LINKER_FLAGS_RELEASE
  ""
)

# TODO: Don't rely on deal.II do stuff ourselves

function(strip_dealii_flags input output)
  string(REGEX REPLACE "-std=[^ ]+" "" stripped "${input}")
  string(STRIP "${stripped}" stripped)
  set(${output} "${stripped}" PARENT_SCOPE)
endfunction()

strip_dealii_flags("${DEAL_II_CXX_FLAGS}" DEAL_II_CXX_FLAGS)
strip_dealii_flags("${DEAL_II_CXX_FLAGS_DEBUG}" DEAL_II_CXX_FLAGS_DEBUG)
strip_dealii_flags("${DEAL_II_CXX_FLAGS_RELEASE}" DEAL_II_CXX_FLAGS_RELEASE)

# TODO: Fix this so it prints with  meaningful variable names
print_variables(
  DEAL_II_CXX_FLAGS
  DEAL_II_CXX_FLAGS_DEBUG
  DEAL_II_CXX_FLAGS_RELEASE
  CMAKE_CXX_FLAGS
  CMAKE_CXX_FLAGS_DEBUG
  CMAKE_CXX_FLAGS_RELEASE
)

separate_arguments(DEAL_II_CXX_FLAGS_LIST UNIX_COMMAND "${DEAL_II_CXX_FLAGS}")
separate_arguments(DEAL_II_CXX_FLAGS_DEBUG_LIST UNIX_COMMAND "${DEAL_II_CXX_FLAGS_DEBUG}")
separate_arguments(
  DEAL_II_CXX_FLAGS_RELEASE_LIST
  UNIX_COMMAND
  "${DEAL_II_CXX_FLAGS_RELEASE}"
)

# Check the deal.II was built with the right dependencies
if(NOT DEAL_II_WITH_MPI)
  message(FATAL_ERROR "deal.II was not built with MPI support!")
endif()
if(NOT DEAL_II_WITH_LAPACK)
  message(FATAL_ERROR "deal.II was not built with LAPACK support!")
endif()
if(NOT DEAL_II_WITH_ZLIB)
  message(FATAL_ERROR "deal.II was not built with zlib support!")
endif()
if(NOT DEAL_II_WITH_P4EST)
  message(FATAL_ERROR "deal.II was not built with p4est support!")
endif()

# Check that deal.II links against the same MPI that we find earlier.
get_target_property(DEAL_II_MPI_LIBS dealii::interface_mpi INTERFACE_LINK_LIBRARIES)
get_target_property(PRISMS_PF_MPI_LIBS MPI::MPI_CXX INTERFACE_LINK_LIBRARIES)

if(NOT DEAL_II_MPI_LIBS STREQUAL PRISMS_PF_MPI_LIBS)
  message(
    WARNING
    "deal.II MPI libraries: ${DEAL_II_MPI_LIBS}\n"
    "Found MPI libraries:   ${PRISMS_PF_MPI_LIBS}\n"
    "These do not match! This may cause runtime issues."
  )
endif()

# Check that deal.II has the same width indices as PRISMS-PF
if(NOT DEAL_II_WITH_64BIT_INDICES STREQUAL PRISMS_PF_64BIT_INDICES)
  message(
    FATAL_ERROR
    "PRISMS-PF and deal.II must both use the same index width. You have:\n"
    "  PRISMS_PF_64BIT_INDICES = ${PRISMS_PF_64BIT_INDICES}\n"
    "  DEAL_II_WITH_64BIT_INDICES = ${DEAL_II_WITH_64BIT_INDICES}\n"
  )
endif()

# TODO: Add caliper
# TODO: Add VTK

# FetchCotent MakeAvailable for the other vendored dependencies
set(PRISMS_PF_VENDORED_PACKAGES libassert)
if(PRISMS_PF_UNIT_TESTS OR PRISMS_PF_REGRESSION_TESTS OR PRISMS_PF_PERFORMANCE_TESTS)
  list(APPEND PRISMS_PF_VENDORED_PACKAGES Catch2)
endif()
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

  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/_deps/libassert-debug/install/include)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/_deps/libassert-release/install/include)

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

  # Still make available other vendored packages (e.g. Catch2)
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
  add_library(prisms_pf_release)

  prisms_pf_configure_target(prisms_pf_debug "Debug")
  prisms_pf_configure_target(prisms_pf_release "Release")

  add_library(prisms_pf::prisms_pf_debug ALIAS prisms_pf_debug)
  add_library(prisms_pf::prisms_pf_release ALIAS prisms_pf_release)

  set(
    PRISMS_PF_TARGETS
    prisms_pf_debug
    prisms_pf_release
  )
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_library(prisms_pf_debug)
  prisms_pf_configure_target(prisms_pf_debug "Debug")
  add_library(prisms_pf::prisms_pf ALIAS prisms_pf_debug)
  add_library(prisms_pf::prisms_pf_debug ALIAS prisms_pf_debug)
  set(PRISMS_PF_TARGETS prisms_pf_debug)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_library(prisms_pf_release)
  prisms_pf_configure_target(prisms_pf_release "Release")
  add_library(prisms_pf::prisms_pf ALIAS prisms_pf_release)
  add_library(prisms_pf::prisms_pf_release ALIAS prisms_pf_release)
  set(PRISMS_PF_TARGETS prisms_pf_release)
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

# =========================================================
# Build exports (find_package with build dir)
# =========================================================

print_subsection("Configuring build-tree export")

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

if(PRISMS_PF_DOCS)
  print_subsection("Configuring PRISMS-PF documentation")

  find_package(Doxygen REQUIRED)

  # Path of Doxyfile and its output
  set(DOXYGEN_INPUT ${CMAKE_SOURCE_DIR}/doc/doxygen/Doxyfile)
  set(DOXYGEN_OUTPUT ${CMAKE_BINARY_DIR}/doc/doxygen/Doxyfile.out)

  # The directory where the documentation will be generated
  set(DOXYGEN_HTML_OUTPUT "${PRISMS_PF_VERSION}")
  set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/doc)

  # Copy input to output
  configure_file(${DOXYGEN_INPUT} ${DOXYGEN_OUTPUT} @ONLY)

  # Add target to generate documentation
  doxygen_add_docs(
    doc
    COMMENT "Generating Doxygen documentation"
    CONFIG_FILE ${DOXYGEN_OUTPUT}
    ALL
  )

  message(
    STATUS
    "Doxygen documentation will be generated using:\n"
    "     ${DOXYGEN_INPUT}\n"
    "   and output to:\n"
    "     ${DOXYGEN_OUTPUT_DIRECTORY}/${DOXYGEN_HTML_OUTPUT}\n"
  )
endif()

# =========================================================
# Final steps to print some nice messages
# =========================================================
