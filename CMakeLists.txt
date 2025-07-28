##
#  CMake for the PRISMS-PF
#  Adapted from the ASPECT CMake file
##

# cmake config
message(STATUS "Using CMake ${CMAKE_VERSION}")
message(STATUS "")
cmake_minimum_required(VERSION 3.13.4)
cmake_policy(VERSION 3.13.4)

# Grab modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

# Load macros
file(GLOB macro_files "${CMAKE_SOURCE_DIR}/cmake/macros/*.cmake")
foreach(file ${macro_files})
  message(STATUS "Include ${file}")
  include(${file})
endforeach()

file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" PRISMS_PF_VERSION LIMIT_COUNT 1)

message(STATUS "")
message(STATUS "=========================================================")
message(STATUS "Configuring PRISMS-PF v${PRISMS_PF_VERSION}")
message(STATUS "=========================================================")
message(STATUS "")

# =========================================================
# Some basic bookkeeping
# =========================================================

# Check that a prior CMakeCache is not located in the build directory
if(EXISTS ${CMAKE_SOURCE_DIR}/CMakeCache.txt)
  message(FATAL_ERROR "Detected the file:\n"
  "${CMAKE_SOURCE_DIR}/CMakeCache.txt\n"
  "in your source directory, which may be leftover from prior builds. " 
  "Please delete the file before running cmake again.")
endif()

# Print the current cmake version
message(STATUS "CMake version: ${CMAKE_VERSION}")

# Grab git information
prisms_pf_git_version()

# Set the standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Setup the build type (debug, release, debugrelease)
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "DebugRelease" 
      CACHE STRING
      "Choose the type of build, options are: Debug, Release and DebugRelease."
      FORCE)
endif()

# Convert build type into the debug and release builds, which may or may 
# not be built.
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR
   "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR
   "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease" )
  message(STATUS "Setting up PRISMS-PF for ${CMAKE_BUILD_TYPE} mode.")
else()
  message(FATAL_ERROR
    "CMAKE_BUILD_TYPE must either be 'Release', 'Debug', or 'DebugRelease', but is set to '${CMAKE_BUILD_TYPE}'.")
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR 
   "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease")
  set(PRISMS_PF_BUILD_DEBUG "ON")
else()
  set(PRISMS_PF_BUILD_DEBUG "OFF")
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR 
   "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease")
  set(PRISMS_PF_BUILD_RELEASE "ON")
else()
  set(PRISMS_PF_BUILD_RELEASE "OFF")
endif()

# Enable autodetection of deal.II dependencies
set(PRISMS_PF_AUTODETECTION ON CACHE BOOL "Autodetection of deal.II dependencies.")
message(STATUS "PRISMS_PF_AUTODETECTION = '${PRISMS_PF_AUTODETECTION}'")

# =========================================================
# External libraries
# =========================================================

message(STATUS "")
message(STATUS "=========================================================")
message(STATUS "Configuring external libraries")
message(STATUS "=========================================================")
message(STATUS "")

# Find deal.II installation
find_package(deal.II 9.6.0 QUIET
	HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR}
  )
if(NOT ${deal.II_FOUND})
  message(FATAL_ERROR "\n*** Could not find a recent version of deal.II. ***\n"
  "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake "
  "or set an environment variable \"DEAL_II_DIR\" that contains a path to a "
  "recent version of deal.II."
  )
endif()

message(STATUS "Found deal.II version ${DEAL_II_PACKAGE_VERSION} at '${deal.II_DIR}'")

set(DEALII_INSTALL_VALID ON)

if(NOT DEAL_II_WITH_LAPACK)
    message(SEND_ERROR
      "\n**deal.II was built without support for lapack!\n"
      )
    set(DEALII_INSTALL_VALID OFF)
endif()

if(NOT DEAL_II_WITH_VTK)
    message(SEND_ERROR
      "\n**deal.II was built without support for VTK!\n"
      )
    set(DEALII_INSTALL_VALID OFF)
endif()


if(NOT DEAL_II_WITH_P4EST)
    message(SEND_ERROR
      "\n**deal.II was built without support for p4est!\n"
      )
    set(DEALII_INSTALL_VALID OFF)
endif()

if(NOT DEAL_II_WITH_MPI)
    message(SEND_ERROR
      "\n**deal.II was built without support for MPI!\n"
      )
    set(DEALII_INSTALL_VALID OFF)
endif()

if(NOT DEALII_INSTALL_VALID)
  message(FATAL_ERROR
    "\nPRISMS-PF requires a deal.II installation with certain features enabled!\n"
    )
endif()

# Load deal.II cached variables
deal_ii_initialize_cached_variables()

# Declare the project
project(prisms_pf CXX)

# Optional deal.II packages

# zlib should be required by p4est in deal.II, so its on by default
set(PRISMS_PF_WITH_ZLIB ON CACHE BOOL "Whether the user wants to compile PRISMS-PF with deal.II's zlib dependency, or not.")
message(STATUS "Using PRISMS_PF_WITH_ZLIB = '${PRISMS_PF_WITH_ZLIB}'")
if(PRISMS_PF_WITH_ZLIB)
  if(DEAL_II_WITH_ZLIB)
    message(STATUS "  Found deal.II installation with zlib")
  else()
    message(FATAL_ERROR "deal.II installation with zlib not found. Disable PRISMS_PF_WITH_ZLIB or recompile deal.II with zlib.")
  endif()
endif()

set(PRISMS_PF_WITH_HDF5 OFF CACHE BOOL "Whether the user wants to compile PRISMS-PF with deal.II's HDF5 dependency, or not.")
if(PRISMS_PF_AUTODETECTION AND DEAL_II_WITH_HDF5)
  set(PRISMS_PF_WITH_HDF5 ON)
endif()
message(STATUS "Using PRISMS_PF_WITH_HDF5 = '${PRISMS_PF_WITH_HDF5}'")
if(PRISMS_PF_WITH_HDF5)
  if(DEAL_II_WITH_HDF5)
    message(STATUS "  Found deal.II installation with HDF5")
  else()
    message(FATAL_ERROR "deal.II installation with HDF5 not found. Disable PRISMS_PF_WITH_HDF5 or recompile deal.II with HDF5.")
  endif()
endif()

set(PRISMS_PF_WITH_CUDA OFF CACHE BOOL "Whether the user wants to compile PRISMS-PF with deal.II's CUDA dependency, or not.")
set(KOKKOS_CUDA_BACKEND OFF CACHE BOOL "Whether the installed Kokkos version has a CUDA backend.")
if(DEAL_II_WITH_KOKKOS)
  find_package(Kokkos QUIET
  HINTS ${KOKKOS_DIR} $ENV{KOKKOS_DIR}
  )
  if(Kokkos_FOUND AND Kokkos_DEVICES MATCHES "CUDA")
    set(KOKKOS_CUDA_BACKEND ON)
  endif()
endif()
if(PRISMS_PF_AUTODETECTION AND KOKKOS_CUDA_BACKEND)
  set(PRISMS_PF_WITH_CUDA ON)
endif()
message(STATUS "Using PRISMS_PF_WITH_CUDA = '${PRISMS_PF_WITH_CUDA}'")
if(PRISMS_PF_WITH_CUDA)
  if(KOKKOS_CUDA_BACKEND)
    message(STATUS "  Found Kokkos installation with CUDA")
  else()
    message(FATAL_ERROR "Kokkos installation with CUDA not found. Disable PRISMS_PF_WITH_CUDA or recompile deal.II with Kokkos that has a CUDA backend.")
  endif()
endif()

set(PRISMS_PF_WITH_SUNDIALS OFF CACHE BOOL "Whether the user wants to compile PRISMS-PF with deal.II's SUNDIALS dependency, or not.")
if(PRISMS_PF_AUTODETECTION AND DEAL_II_WITH_SUNDIALS)
  set(PRISMS_PF_WITH_SUNDIALS ON)
endif()
message(STATUS "Using PRISMS_PF_WITH_SUNDIALS = '${PRISMS_PF_WITH_SUNDIALS}'")
if(PRISMS_PF_WITH_SUNDIALS)
  if(DEAL_II_WITH_SUNDIALS)
    message(STATUS "  Found deal.II installation with SUNDIALS")
  else()
    message(FATAL_ERROR "deal.II installation with SUNDIALS not found. Disable PRISMS_PF_WITH_SUNDIALS or recompile deal.II with SUNDIALS.")
  endif()
endif()

# Caliper
set(PRISMS_PF_WITH_CALIPER OFF CACHE BOOL "Whether the user wants to compile PRISMS-PF with the profiling code Caliper, or not.")
if(PRISMS_PF_AUTODETECTION)
  find_package(CALIPER QUIET
  HINTS ${CALIPER_DIR} $ENV{CALIPER_DIR}
  )
  if(${CALIPER_FOUND})
    set(PRISMS_PF_WITH_CALIPER ON)
  endif()
endif()
message(STATUS "Using PRISMS_PF_WITH_CALIPER = '${PRISMS_PF_WITH_CALIPER}'")
if(PRISMS_PF_WITH_CALIPER)
  find_package(CALIPER QUIET
    HINTS ${CALIPER_DIR} $ENV{CALIPER_DIR}
    )
  message(STATUS "CALIPER_DIR: ${CALIPER_DIR}")
  if(${CALIPER_FOUND})
    include_directories(${CALIPER_INCLUDE_DIR})
    message(STATUS "  Caliper found at ${CALIPER_DIR}")
  else()
    message(FATAL_ERROR "Caliper not found. Disable PRISMS_PF_WITH_CALIPER or specify a hint to your installation directory with CALIPER_DIR")
  endif()
endif()

# Create compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(FORCE_COLORED_OUTPUT ON CACHE BOOL "Forces colored output when compiling with gcc and clang.")

# =========================================================
# Configure PRISMS-PF Targets
# =========================================================

message(STATUS "")
message(STATUS "=========================================================")
message(STATUS "Configuring PRISMS-PF build targets")
message(STATUS "=========================================================")
message(STATUS "")

# Grab flag for additional optimizations
set(ADDITIONAL_OPTIMIZATIONS OFF CACHE BOOL "Whether the user wants to enable additional optimizations, or not.")
message(STATUS "ADDITIONAL_OPTIMIZATIONS = '${ADDITIONAL_OPTIMIZATIONS}'")
if(ADDITIONAL_OPTIMIZATIONS)
  message(WARNING "ADDITIONAL_OPTIMIZATIONS is enabled. Procede with caution.")
endif()

# Generate config.h to enable and disable certain features within the source code.
set(PRISMS_PF_SOURCE_DIR ${CMAKE_SOURCE_DIR})

configure_file(
  ${CMAKE_SOURCE_DIR}/include/prismspf/config.h.in
  ${CMAKE_BINARY_DIR}/include/prismspf/config.h
  )

# Collect source files
add_subdirectory(src)

# Test stuff goes here

# Make and ninja build options
if(CMAKE_GENERATOR MATCHES "Ninja")
  set(_make_command "$ ninja")
else()
  set(_make_command "$ make")
endif()

# Debug and release targets
if(${DEAL_II_BUILD_TYPE} MATCHES "DebugRelease")
  add_custom_target(release
    COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Release .
    COMMAND ${CMAKE_COMMAND} -E echo "***"
    COMMAND ${CMAKE_COMMAND} -E echo "*** Switched to Release mode. Now recompile with: ${_make_command}"
    COMMAND ${CMAKE_COMMAND} -E echo "***"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
    COMMENT "switching to RELEASE mode..."
  )

  add_custom_target(debug
    COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Debug .
    COMMAND ${CMAKE_COMMAND} -E echo "***"
    COMMAND ${CMAKE_COMMAND} -E echo "*** Switched to Debug mode. Now recompile with: ${_make_command}"
    COMMAND ${CMAKE_COMMAND} -E echo "***"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
    COMMENT "switching to DEBUG mode..."
  )

  add_custom_target(debugrelease
    COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=DebugRelease .
    COMMAND ${CMAKE_COMMAND} -E echo "***"
    COMMAND ${CMAKE_COMMAND} -E echo "*** Switched to Debug and Release mode. Now recompile with: ${_make_command}"
    COMMAND ${CMAKE_COMMAND} -E echo "***"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
    COMMENT "switching to DEBUG/RELEASE mode..."
  )
endif()

# Add distclean target to clean build
add_custom_target(distclean
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
  COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
  COMMAND ${CMAKE_COMMAND} -E remove
    CMakeCache.txt cmake_install.cmake Makefile
    build.ninja rules.ninja .ninja_deps .ninja_log
    clang-tidy.log cppcheck.log cppcheck_unused.log output.txt error.txt
    summary.log detailed.log
  COMMENT "distclean invoked"
  )


file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
"message(
\"###
#
#  PRISMS-PF set up with ${DEAL_II_PACKAGE_NAME}-${DEAL_II_PACKAGE_VERSION} found at
#      ${DEAL_II_PATH}
#
#  CMAKE_BUILD_TYPE:          ${CMAKE_BUILD_TYPE}
#
#  You can now run
#      ${_make_command}                - to compile and link
#      ${_make_command} debug          - to switch the build type to 'Debug'
#      ${_make_command} release        - to switch the build type to 'Release'
#      ${_make_command} debugrelease   - to switch the build type to compile both
#      ${_make_command} clean          - to remove the generated executable as well as
#                               all intermediate compilation files
#      ${_make_command} distclean      - to clean the directory from all generated
#                               files (includes clean, runclean and the removal
#                               of the generated build system)
#      ${_make_command} info           - to view this message again
\")")

add_custom_target(info
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
  )

if (${FORCE_COLORED_OUTPUT})
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER MATCHES "AppleClang")
    string(APPEND DEAL_II_CXX_FLAGS_DEBUG " -fcolor-diagnostics")
    string(APPEND DEAL_II_CXX_FLAGS_RELEASE " -fcolor-diagnostics")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    string(APPEND DEAL_II_CXX_FLAGS_DEBUG " -fdiagnostics-color=always")
    string(APPEND DEAL_II_CXX_FLAGS_RELEASE " -fdiagnostics-color=always")
  endif()
endif()

# deal.II versions >=9.5 disable deprecation warnings in user code. Enable
# the warnings again by removing the flag that disables them.
string(REPLACE "-Wno-deprecated-declarations" "" DEAL_II_WARNING_FLAGS "${DEAL_II_WARNING_FLAGS}")

# Set additional compiler flags
set(PRISMS_PF_ADDITIONAL_CXX_FLAGS "" CACHE STRING "Additional CMAKE_CXX_FLAGS applied after the deal.II options.")

if(NOT PRISMS_PF_ADDITIONAL_CXX_FLAGS STREQUAL "")
  message(STATUS "Appending PRISMS_PF_ADDITIONAL_CXX_FLAGS: '${PRISMS_PF_ADDITIONAL_CXX_FLAGS}':")
  string(APPEND DEAL_II_CXX_FLAGS_DEBUG " ${PRISMS_PF_ADDITIONAL_CXX_FLAGS}")
  string(APPEND DEAL_II_CXX_FLAGS_RELEASE " ${PRISMS_PF_ADDITIONAL_CXX_FLAGS}")
  message(STATUS "  DEAL_II_WARNING_FLAGS: ${DEAL_II_WARNING_FLAGS}")
  message(STATUS "  DEAL_II_CXX_FLAGS_DEBUG: ${DEAL_II_CXX_FLAGS_DEBUG}")
  message(STATUS "  DEAL_II_CXX_FLAGS_RELEASE: ${DEAL_II_CXX_FLAGS_RELEASE}")
endif()

if(${PRISMS_PF_BUILD_DEBUG} STREQUAL "ON")
  set(LIBRARY_NAME_DEBUG "prisms_pf_debug")
  add_library(${LIBRARY_NAME_DEBUG} STATIC ${PRISMS_PF_SOURCE_FILES})
  set_property(TARGET ${LIBRARY_NAME_DEBUG} PROPERTY OUTPUT_NAME prisms-pf-debug)
  target_include_directories(${LIBRARY_NAME_DEBUG} PRIVATE ${CMAKE_BINARY_DIR}/include include)
  deal_ii_setup_target(${LIBRARY_NAME_DEBUG} DEBUG)

  if(${PRISMS_PF_WITH_CALIPER})
    target_link_libraries(${LIBRARY_NAME_DEBUG} caliper)
  endif()
endif()

if(${PRISMS_PF_BUILD_RELEASE} STREQUAL "ON")
  set(LIBRARY_NAME_RELEASE "prisms_pf_release")
  add_library(${LIBRARY_NAME_RELEASE} STATIC ${PRISMS_PF_SOURCE_FILES})
  set_property(TARGET ${LIBRARY_NAME_RELEASE} PROPERTY OUTPUT_NAME prisms-pf-release)
  target_include_directories(${LIBRARY_NAME_RELEASE} PRIVATE ${CMAKE_BINARY_DIR}/include include)
  deal_ii_setup_target(${LIBRARY_NAME_RELEASE} RELEASE)

  if(${PRISMS_PF_WITH_CALIPER})
    target_link_libraries(${LIBRARY_NAME_RELEASE} caliper)
  endif()
endif()

# Unwrap the compiler compile_commands.json after the build
set(UNWRAP_COMPILER OFF CACHE BOOL "Whether the user wants to unwrap the compiler, or not.")
if(${UNWRAP_COMPILER} STREQUAL "ON")
  add_custom_command(
    TARGET ${LIBRARY_NAME_DEBUG} POST_BUILD
    COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/contrib/utilities/unwrap_compile_commands.sh
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Unwrapping compiler for compile_commands.json"
  )
endif()

# Write config to file
include(${CMAKE_SOURCE_DIR}/cmake/write_config.cmake)

# Generate configuration file for applications
configure_file(
  ${CMAKE_SOURCE_DIR}/cmake/prisms_pf_config.cmake.in
  ${CMAKE_BINARY_DIR}/cmake/prisms_pf_config.cmake
  @ONLY
)

# Print the config
file(READ ${CMAKE_SOURCE_DIR}/summary.log PRISMS_PF_SUMMARY)
message("${PRISMS_PF_SUMMARY}")
