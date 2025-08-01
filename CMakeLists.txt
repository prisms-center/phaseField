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

# Check that a prior CMakeCache is not located in the build directory
if(EXISTS ${CMAKE_SOURCE_DIR}/CMakeCache.txt)
    message(
        FATAL_ERROR
        "Detected the file:\n"
        "${CMAKE_SOURCE_DIR}/CMakeCache.txt\n"
        "in your source directory, which may be leftover from prior builds. "
        "Please delete the file before running cmake again."
    )
endif()

# Grab git information
prisms_pf_git_version()

# Include configuration options
include(${CMAKE_SOURCE_DIR}/cmake/setup_cached_variables.cmake)

# Declare the project
project(prisms_pf CXX)

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
set(PRISMS_PF_MODULES DEAL_II VTK CALIPER)
foreach(_module ${PRISMS_PF_MODULES})
    include(${CMAKE_SOURCE_DIR}/cmake/modules/FindPRISMS_PF_${_module}.cmake)
endforeach()

message(STATUS "Using PRISMS_PF_WITH_ZLIB = '${PRISMS_PF_WITH_ZLIB}'")
message(STATUS "Using PRISMS_PF_WITH_HDF5 = '${PRISMS_PF_WITH_HDF5}'")
message(STATUS "Using PRISMS_PF_WITH_SUNDIALS = '${PRISMS_PF_WITH_SUNDIALS}'")
message(STATUS "Using PRISMS_PF_WITH_CUDA = '${PRISMS_PF_WITH_CUDA}'")
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
set(PRISMS_PF_SOURCE_DIR ${CMAKE_SOURCE_DIR})

if(ADDITIONAL_DEGREES)
    set(EXPAND_DEGREES "; 4; 5; 6")
endif()

configure_file(
    ${CMAKE_SOURCE_DIR}/include/prismspf/config.h.in
    ${CMAKE_BINARY_DIR}/include/prismspf/config.h
)

# Configure the template arguments file
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/templates.in
    ${CMAKE_BINARY_DIR}/cmake/templates
)

# Add the script files
add_subdirectory("${CMAKE_SOURCE_DIR}/cmake/scripts")

# Collect source files
add_subdirectory(src)

# Grab the source and include files
get_property(PRISMS_PF_SOURCE_FILES GLOBAL PROPERTY PRISMS_PF_SOURCE_FILES)
get_property(PRISMS_PF_HEADER_FILES GLOBAL PROPERTY PRISMS_PF_HEADER_FILES)

# Test stuff goes here

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
        COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Release .
        COMMAND ${CMAKE_COMMAND} -E echo "***"
        COMMAND
            ${CMAKE_COMMAND} -E echo
            "*** Switched to Release mode. Now recompile with: ${_make_command}"
        COMMAND ${CMAKE_COMMAND} -E echo "***"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        VERBATIM
        COMMENT "switching to RELEASE mode..."
    )

    add_custom_target(
        debug
        COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Debug .
        COMMAND ${CMAKE_COMMAND} -E echo "***"
        COMMAND
            ${CMAKE_COMMAND} -E echo
            "*** Switched to Debug mode. Now recompile with: ${_make_command}"
        COMMAND ${CMAKE_COMMAND} -E echo "***"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        VERBATIM
        COMMENT "switching to DEBUG mode..."
    )

    add_custom_target(
        debugrelease
        COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=DebugRelease .
        COMMAND ${CMAKE_COMMAND} -E echo "***"
        COMMAND
            ${CMAKE_COMMAND} -E echo
            "*** Switched to Debug and Release mode. Now recompile with: ${_make_command}"
        COMMAND ${CMAKE_COMMAND} -E echo "***"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        VERBATIM
        COMMENT "switching to DEBUG/RELEASE mode..."
    )
endif()

# Add distclean target to clean build
add_custom_target(
    distclean
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
    COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
    COMMAND
        ${CMAKE_COMMAND} -E remove CMakeCache.txt cmake_install.cmake Makefile
        build.ninja rules.ninja .ninja_deps .ninja_log clang-tidy.log
        cppcheck.log cppcheck_unused.log output.txt error.txt summary.log
        detailed.log
    COMMENT "distclean invoked"
)

file(
    WRITE
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
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
        ${CMAKE_COMMAND} -P
        ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
)

if(${FORCE_COLORED_OUTPUT})
    if(
        CMAKE_CXX_COMPILER_ID MATCHES "Clang"
        OR CMAKE_CXX_COMPILER MATCHES "AppleClang"
    )
        string(APPEND DEAL_II_CXX_FLAGS_DEBUG " -fcolor-diagnostics")
        string(APPEND DEAL_II_CXX_FLAGS_RELEASE " -fcolor-diagnostics")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        string(APPEND DEAL_II_CXX_FLAGS_DEBUG " -fdiagnostics-color=always")
        string(APPEND DEAL_II_CXX_FLAGS_RELEASE " -fdiagnostics-color=always")
    endif()
endif()

# deal.II versions >=9.5 disable deprecation warnings in user code. Enable
# the warnings again by removing the flag that disables them.
string(
    REPLACE
    "-Wno-deprecated-declarations"
    ""
    DEAL_II_WARNING_FLAGS
    "${DEAL_II_WARNING_FLAGS}"
)

# Set additional compiler flags
set(PRISMS_PF_ADDITIONAL_CXX_FLAGS
    ""
    CACHE STRING
    "Additional CMAKE_CXX_FLAGS applied after the deal.II options."
)

if(NOT PRISMS_PF_ADDITIONAL_CXX_FLAGS STREQUAL "")
    message(
        STATUS
        "Appending PRISMS_PF_ADDITIONAL_CXX_FLAGS: '${PRISMS_PF_ADDITIONAL_CXX_FLAGS}':"
    )
    string(APPEND DEAL_II_CXX_FLAGS_DEBUG " ${PRISMS_PF_ADDITIONAL_CXX_FLAGS}")
    string(
        APPEND
        DEAL_II_CXX_FLAGS_RELEASE
        " ${PRISMS_PF_ADDITIONAL_CXX_FLAGS}"
    )
    message(STATUS "  DEAL_II_WARNING_FLAGS: ${DEAL_II_WARNING_FLAGS}")
    message(STATUS "  DEAL_II_CXX_FLAGS_DEBUG: ${DEAL_II_CXX_FLAGS_DEBUG}")
    message(STATUS "  DEAL_II_CXX_FLAGS_RELEASE: ${DEAL_II_CXX_FLAGS_RELEASE}")
endif()

if(${PRISMS_PF_BUILD_DEBUG} STREQUAL "ON")
    set(LIBRARY_NAME_DEBUG "prisms_pf_debug")
    add_library(${LIBRARY_NAME_DEBUG} STATIC)
    target_sources(${LIBRARY_NAME_DEBUG} PRIVATE ${PRISMS_PF_SOURCE_FILES})
    set_property(
        TARGET ${LIBRARY_NAME_DEBUG}
        PROPERTY OUTPUT_NAME prisms-pf-debug
    )
    target_link_libraries(
        ${LIBRARY_NAME_DEBUG}
        prisms_pf_core_debug
        prisms_pf_solvers_debug
        prisms_pf_user_inputs_debug
        prisms_pf_utilities_debug
    )
    # Add explicit dependencies to ensure proper build order with Ninja
    add_dependencies(
        ${LIBRARY_NAME_DEBUG}
        prisms_pf_core_debug
        prisms_pf_solvers_debug
        prisms_pf_user_inputs_debug
        prisms_pf_utilities_debug
    )
    target_include_directories(
        ${LIBRARY_NAME_DEBUG}
        PRIVATE ${CMAKE_BINARY_DIR}/include ${CMAKE_BINARY_DIR}/src include
    )

    if(${VTK_BUILT_SEPARATELY})
        include_directories(SYSTEM ${VTK_INCLUDE_DIR})
        target_link_libraries(${LIBRARY_NAME_DEBUG} ${VTK_NEW_LIBRARIES})
    endif()
    if(${PRISMS_PF_WITH_CALIPER})
        target_link_libraries(${LIBRARY_NAME_DEBUG} caliper)
    endif()

    deal_ii_setup_target(${LIBRARY_NAME_DEBUG} DEBUG)
endif()

if(${PRISMS_PF_BUILD_RELEASE} STREQUAL "ON")
    set(LIBRARY_NAME_RELEASE "prisms_pf_release")
    add_library(${LIBRARY_NAME_RELEASE} STATIC)
    target_sources(${LIBRARY_NAME_RELEASE} PRIVATE ${PRISMS_PF_SOURCE_FILES})
    set_property(
        TARGET ${LIBRARY_NAME_RELEASE}
        PROPERTY OUTPUT_NAME prisms-pf-release
    )
    target_link_libraries(
        ${LIBRARY_NAME_RELEASE}
        prisms_pf_core_release
        prisms_pf_solvers_release
        prisms_pf_user_inputs_release
        prisms_pf_utilities_release
    )
    # Add explicit dependencies to ensure proper build order with Ninja
    add_dependencies(
        ${LIBRARY_NAME_RELEASE}
        prisms_pf_core_release
        prisms_pf_solvers_release
        prisms_pf_user_inputs_release
        prisms_pf_utilities_release
    )
    target_include_directories(
        ${LIBRARY_NAME_RELEASE}
        PRIVATE ${CMAKE_BINARY_DIR}/include ${CMAKE_BINARY_DIR}/src include
    )

    if(${VTK_BUILT_SEPARATELY})
        include_directories(SYSTEM ${VTK_INCLUDE_DIR})
        target_link_libraries(${LIBRARY_NAME_RELEASE} ${VTK_NEW_LIBRARIES})
    endif()
    if(${PRISMS_PF_WITH_CALIPER})
        target_link_libraries(${LIBRARY_NAME_RELEASE} caliper)
    endif()

    deal_ii_setup_target(${LIBRARY_NAME_RELEASE} RELEASE)
endif()

# Unwrap the compiler compile_commands.json after the build
if(${UNWRAP_COMPILER} STREQUAL "ON")
    add_custom_command(
        TARGET ${LIBRARY_NAME_DEBUG}
        POST_BUILD
        COMMAND
            /bin/bash
            ${CMAKE_SOURCE_DIR}/contrib/utilities/unwrap_compile_commands.sh
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
