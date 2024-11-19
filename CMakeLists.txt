##
#  CMake for the PRISMS-PF
#  Adapted from the ASPECT CMake file
##

cmake_minimum_required(VERSION 3.1.0)

project(prisms_pf)

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

# =========================================================
# External libraries
# =========================================================

message(STATUS "")
message(STATUS "=========================================================")
message(STATUS "Configuring external libraries")
message(STATUS "=========================================================")
message(STATUS "")

# Find deal.II installation
find_package(deal.II 9.2.0 QUIET
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

if(NOT DEAL_II_WITH_P4EST)
    message(SEND_ERROR
      "\n**deal.II was built without support for p4est!\n"
      )
    set(DEALII_INSTALL_VALID OFF)
endif()

if(NOT DEALII_INSTALL_VALID)
  message(FATAL_ERROR
    "\nPRISMS-PF requires a deal.II installation with certain features enabled!\n"
    )
endif()

deal_ii_initialize_cached_variables()

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

# Generate config.h to enable and disable certain features within the source code.
set(PRISMS_PF_SOURCE_DIR ${CMAKE_SOURCE_DIR})

configure_file(
  ${CMAKE_SOURCE_DIR}/include/config.h.in
  ${CMAKE_BINARY_DIR}/include/config.h
  )

# Collect source files
file(GLOB_RECURSE PRISMS_PF_SOURCE_FILES "src/*.cc" "include/*.h")
file(GLOB_RECURSE PRISMS_PF_UNIT_TEST_FILES "unit_tests/*.cc" "contrib/catch/catch.hpp")

# Setup include directories
include_directories(BEFORE ${CMAKE_BINARY_DIR}/include include)
include_directories(AFTER contrib/catch)

# Test stuff goes here

# Make and ninja build options
if(CMAKE_GENERATOR MATCHES "Ninja")
  set(_make_command "$ ninja")
else()
  set(_make_command " $ make")
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
  COMMENT "distclean invoked"
  )


file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
"message(
\"###
#
#  Project ${TARGET_EXE} set up with  ${DEAL_II_PACKAGE_NAME}-${DEAL_II_PACKAGE_VERSION}  found at
#      ${DEAL_II_PATH}
#
#  CMAKE_BUILD_TYPE:          ${CMAKE_BUILD_TYPE}
#
#  You can now run
#      ${_make_command}                - to compile and link ${TARGET_EXE}
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
endif()

if(${PRISMS_PF_BUILD_RELEASE} STREQUAL "ON")
  set(LIBRARY_NAME_RELEASE "prisms_pf_release")
  add_library(${LIBRARY_NAME_RELEASE} STATIC ${PRISMS_PF_SOURCE_FILES})
  set_property(TARGET ${LIBRARY_NAME_RELEASE} PROPERTY OUTPUT_NAME prisms-pf-release)
  target_include_directories(${LIBRARY_NAME_RELEASE} PRIVATE ${CMAKE_BINARY_DIR}/include include)
  deal_ii_setup_target(${LIBRARY_NAME_RELEASE} RELEASE)
endif()