##
#  CMake for the PRISMS-PF core library
##

# cmake config
message(STATUS "Using CMake ${CMAKE_VERSION}")
message(STATUS "")
cmake_minimum_required(VERSION 3.18)
cmake_policy(VERSION 3.18)

# Grab modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

# Load macros
file(GLOB macro_files "${CMAKE_SOURCE_DIR}/cmake/macros/*.cmake")
foreach(file ${macro_files})
  message(STATUS "Include ${file}")
  include(${file})
endforeach()

# Grab the version of PRISMS-PF
file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" PRISMS_PF_VERSION LIMIT_COUNT 1)

message(STATUS "")
message(STATUS "=========================================================")
message(STATUS "Configuring PRISMS-PF v${PRISMS_PF_VERSION}")
message(STATUS "=========================================================")
message(STATUS "")

# =========================================================
# Some basic bookkeeping
# =========================================================

# Grab git information
prisms_pf_git_version()

# Include configuration options
include(${CMAKE_SOURCE_DIR}/cmake/setup_cached_variables.cmake)

# Declare the project
project(prisms_pf CXX)

# Set the default install prefix if it was not set by the user
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(
    CMAKE_INSTALL_PREFIX
    "${CMAKE_BINARY_DIR}/install"
    CACHE PATH
    "Installation prefix"
    FORCE
  )
endif()

# Include the PRISMS-PF setup
include(${CMAKE_SOURCE_DIR}/cmake/setup_prisms_pf.cmake)

# =========================================================
# External libraries
# =========================================================

message(STATUS "")
message(STATUS "=========================================================")
message(STATUS "Configuring external libraries")
message(STATUS "=========================================================")
message(STATUS "")

# Process each of the modules
set(
  PRISMS_PF_MODULES
  DEAL_II
  VTK
  CALIPER
)
foreach(_module ${PRISMS_PF_MODULES})
  include(${CMAKE_SOURCE_DIR}/cmake/modules/FindPRISMS_PF_${_module}.cmake)
endforeach()

message(STATUS "Using PRISMS_PF_WITH_ZLIB = '${PRISMS_PF_WITH_ZLIB}'")
message(STATUS "Using PRISMS_PF_WITH_VTK = '${PRISMS_PF_WITH_VTK}'")
message(STATUS "Using PRISMS_PF_WITH_CALIPER = '${PRISMS_PF_WITH_CALIPER}'")

# =========================================================
# Configure PRISMS-PF Targets
# =========================================================

message(STATUS "")
message(STATUS "=========================================================")
message(STATUS "Configuring PRISMS-PF build targets")
message(STATUS "=========================================================")
message(STATUS "")

# Generate config.h to enable and disable certain features within the source code.
set(PRISMS_PF_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}")

if(ADDITIONAL_DEGREES)
  set(EXPAND_DEGREES "; 3; 4; 5; 6")
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/prismspf/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/prismspf/config.h
)

# Configure the template arguments file
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/templates.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/templates
)

# Add the script files
add_subdirectory("${CMAKE_SOURCE_DIR}/cmake/scripts")

# Collect source files
add_subdirectory(src)

# Grab the source and include files
get_property(PRISMS_PF_SOURCE_FILES GLOBAL PROPERTY PRISMS_PF_SOURCE_FILES)
get_property(PRISMS_PF_HEADER_FILES GLOBAL PROPERTY PRISMS_PF_HEADER_FILES)

# Make and ninja build options
if(CMAKE_GENERATOR MATCHES "Ninja")
  set(_make_command "$ ninja")
else()
  set(_make_command "$ make")
endif()

# Debug and release targets
if(${DEAL_II_BUILD_TYPE} MATCHES "DebugRelease")
  add_custom_target(
    release
    COMMAND
      ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Release .
    COMMAND
      ${CMAKE_COMMAND} -E echo "***"
    COMMAND
      ${CMAKE_COMMAND} -E echo
      "*** Switched to Release mode. Now recompile with: ${_make_command}"
    COMMAND
      ${CMAKE_COMMAND} -E echo "***"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
    COMMENT "switching to RELEASE mode..."
  )

  add_custom_target(
    debug
    COMMAND
      ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Debug .
    COMMAND
      ${CMAKE_COMMAND} -E echo "***"
    COMMAND
      ${CMAKE_COMMAND} -E echo
      "*** Switched to Debug mode. Now recompile with: ${_make_command}"
    COMMAND
      ${CMAKE_COMMAND} -E echo "***"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
    COMMENT "switching to DEBUG mode..."
  )

  add_custom_target(
    debugrelease
    COMMAND
      ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=DebugRelease .
    COMMAND
      ${CMAKE_COMMAND} -E echo "***"
    COMMAND
      ${CMAKE_COMMAND} -E echo
      "*** Switched to Debug and Release mode. Now recompile with: ${_make_command}"
    COMMAND
      ${CMAKE_COMMAND} -E echo "***"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
    COMMENT "switching to DEBUG/RELEASE mode..."
  )
endif()

# Add distclean target to clean build
add_custom_target(
  distclean
  COMMAND
    ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
  COMMAND
    ${CMAKE_COMMAND} -E remove_directory CMakeFiles
  COMMAND
    ${CMAKE_COMMAND} -E remove CMakeCache.txt cmake_install.cmake Makefile build.ninja
    rules.ninja .ninja_deps .ninja_log clang-tidy.log cppcheck.log cppcheck_unused.log
    output.txt error.txt summary.log detailed.log
  COMMENT "distclean invoked"
)

file(
  WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
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
\")"
)

add_custom_target(
  info
  COMMAND
    ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
)

# Set additional compiler flags
set(
  PRISMS_PF_ADDITIONAL_CXX_FLAGS
  ""
  CACHE STRING
  "Additional CMAKE_CXX_FLAGS applied after the deal.II options."
)

if(NOT PRISMS_PF_ADDITIONAL_CXX_FLAGS STREQUAL "")
  message(
    STATUS
    "Appending PRISMS_PF_ADDITIONAL_CXX_FLAGS: '${PRISMS_PF_ADDITIONAL_CXX_FLAGS}':"
  )
  append_flags(PRISMS_PF_ADDITIONAL_CXX_FLAGS PRISMS_PF_CXX_FLAGS)
  message(STATUS "  PRISMS_PF_WARNING_FLAGS: ${PRISMS_PF_WARNING_FLAGS}")
  message(STATUS "  PRISMS_PF_CXX_FLAGS: ${PRISMS_PF_CXX_FLAGS}")
  message(STATUS "  PRISMS_PF_CXX_FLAGS_DEBUG: ${PRISMS_PF_CXX_FLAGS_DEBUG}")
  message(STATUS "  PRISMS_PF_CXX_FLAGS_RELEASE: ${PRISMS_PF_CXX_FLAGS_RELEASE}")
endif()

foreach(_build ${PRISMS_PF_BUILD_TYPES})
  # Create a target for each build type. Note that the build types are
  # mixed case (e.g., Debug) so we need to convert them to lowercase and
  # uppercase for various uses.
  string(TOLOWER ${_build} _build_lowercase)
  string(TOUPPER ${_build} _build_uppercase)
  set(_library "prisms_pf_${_build_lowercase}")
  list(APPEND PRISMS_PF_TARGETS ${_library})

  # Add the library and set some properties
  add_library(${_library} STATIC)
  set_target_properties(
    ${_target}
    PROPERTIES
      LINKER_LANGUAGE
        CXX
      OUTPUT_NAME
        "prisms-pf-${_build_lowercase}"
  )
  target_sources(${_library} PRIVATE ${PRISMS_PF_SOURCE_FILES})

  # Collect the libraries to link against and the include directories
  set(
    _library_dependencies
    prisms_pf_core_${_build_lowercase}
    prisms_pf_solvers_${_build_lowercase}
    prisms_pf_user_inputs_${_build_lowercase}
    prisms_pf_utilities_${_build_lowercase}
  )
  list(APPEND PRISMS_PF_TARGETS ${_library_dependencies})

  # Link the PRISMS-PF libraries
  target_link_libraries(${_library} PUBLIC ${_library_dependencies})
  # Add explicit dependencies to ensure proper build order with Ninja
  add_dependencies(${_library} ${_library_dependencies})

  # Set the include directories
  target_include_directories(
    ${_library}
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:include>
    PRIVATE ${CMAKE_SOURCE_DIR}/src
  )

  target_compile_options(
    ${_library}
    PRIVATE
      $<$<COMPILE_LANGUAGE:CXX>:
      ${PRISMS_PF_WARNING_FLAGS}
      ${PRISMS_PF_CXX_FLAGS}
      ${PRISMS_PF_CXX_FLAGS_${_build_uppercase}}>
  )
  target_link_options(
    ${_library}
    PRIVATE
      $<$<LINK_LANGUAGE:CXX>:
      ${PRISMS_PF_LINKER_FLAGS}
      ${PRISMS_PF_LINKER_FLAGS_${_build_uppercase}}>
  )

  # Add other dependencies, making sure they are public so that they
  # propagate to targets that link against this library

  # VTK
  if(${VTK_BUILT_SEPARATELY})
    target_include_directories(${_library} SYSTEM PUBLIC ${VTK_INCLUDE_DIR})
    target_link_libraries(${_library} PUBLIC ${VTK_NEW_LIBRARIES})
  endif()

  # caliper
  if(${PRISMS_PF_WITH_CALIPER})
    target_link_libraries(${_library} PUBLIC caliper)
  endif()

  # deal.II
  target_link_libraries(${_library} PUBLIC ${DEAL_II_TARGET_${_build_uppercase}})
endforeach()

# =========================================================
# Installation
# =========================================================
include(GNUInstallDirs)

# Install the libraries
install(
  TARGETS
    ${PRISMS_PF_TARGETS}
  EXPORT prisms_pf-targets
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install the headers
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
  DESTINATION include
  FILES_MATCHING
  PATTERN "*.h"
  PATTERN "*.hpp"
)

# Install the generated config.h
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/include/prismspf/config.h
    DESTINATION include/prismspf
)

# Install the cmake folder
install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cmake/
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prisms_pf
    FILES_MATCHING
    PATTERN "*.cmake"
    PATTERN "*.cc"
    PATTERN "CMakeLists.txt"
    PATTERN "*.cmake.in" EXCLUDE
)

# Install the generated Config file
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/prisms_pf-config.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prisms_pf
)

# Install the CMake config files
install(
  EXPORT prisms_pf-targets
  FILE prisms_pf-targets.cmake
  NAMESPACE prisms_pf::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prisms_pf
)

# Install the VERSION file
install(
  FILES ${CMAKE_CURRENT_SOURCE_DIR}/VERSION
  DESTINATION .
)

# Install the LICENSE file
install(
  FILES ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE
  DESTINATION .
)

# Install the cmake/template file
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/cmake/templates"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prisms_pf
)

# =========================================================
# Testing
# =========================================================

# Unit tests
if(UNIT_TESTS)
  enable_testing()
  add_subdirectory(tests/unit_tests)

  add_custom_target(
    run_tests
    COMMAND
      ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS
      test_core
      test_field_input
      test_user_inputs
      test_utilities
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Run all tests"
  )
endif()

# Regression tests
if(REGRESSION_TESTS)
  # Compile each of the regression tests and run the python script
  add_subdirectory(tests/regression_tests)
endif()

# Applications
if(BUILD_APPLICATIONS)
  # Compile each of the applications
  add_subdirectory(applications)
endif()

# =========================================================
# Final steps
# =========================================================

# Write config to file
include(${CMAKE_SOURCE_DIR}/cmake/write_config.cmake)

# Generate configuration file for applications
include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/prisms_pf-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/prisms_pf-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prisms_pf
)

# Unwrap the compiler compile_commands.json after the build
if(${UNWRAP_COMPILER} STREQUAL "ON")
  add_custom_command(
    TARGET prisms_pf_debug
    POST_BUILD
    COMMAND
      /bin/bash ${CMAKE_SOURCE_DIR}/contrib/utilities/unwrap_compile_commands.sh
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Unwrapping compiler for compile_commands.json"
  )
endif()

# Export targets
list(REMOVE_DUPLICATES PRISMS_PF_TARGETS)
export(
  TARGETS
    ${PRISMS_PF_TARGETS}
  FILE "${CMAKE_BINARY_DIR}/cmake/prisms_pf-targets.cmake"
  NAMESPACE prisms_pf::
)

# Print the config
file(READ ${CMAKE_SOURCE_DIR}/summary.log PRISMS_PF_SUMMARY)
message("${PRISMS_PF_SUMMARY}")
