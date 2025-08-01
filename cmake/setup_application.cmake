# General setup for a PRISMS-PF application

# Check that PRISMS_PF_CORE_DIR is set
if(PRISMS_PF_CORE_DIR STREQUAL "")
    message(FATAL_ERROR "PRISMS_PF_CORE_DIR is not set")
endif()

# Load macros
file(GLOB macro_files "${PRISMS_PF_CORE_DIR}/cmake/macros/*.cmake")
foreach(file ${macro_files})
    message(STATUS "Include ${file}")
    include(${file})
endforeach()

# Grab the version of PRISMS-PF
file(STRINGS "${PRISMS_PF_CORE_DIR}/VERSION" PRISMS_PF_VERSION LIMIT_COUNT 1)

message(STATUS "")
message(STATUS "=========================================================")
message(STATUS "Configuring application with PRISMS-PF v${PRISMS_PF_VERSION}")
message(STATUS "=========================================================")
message(STATUS "")

# Set the standard to C++20
set_cpp_standard(20)

# Setup the build type (debug, release, debugrelease)
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE
        "DebugRelease"
        CACHE STRING
        "Choose the type of build, options are: Debug, Release and DebugRelease."
        FORCE
    )
endif()

# Convert build type into the debug and release builds, which may or may
# not be built.
if(
    "${CMAKE_BUILD_TYPE}" STREQUAL "Release"
    OR "${CMAKE_BUILD_TYPE}" STREQUAL "Debug"
    OR "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease"
)
    message(
        STATUS
        "Setting up PRISMS-PF application for ${CMAKE_BUILD_TYPE} mode."
    )
else()
    message(
        FATAL_ERROR
        "CMAKE_BUILD_TYPE must either be 'Release', 'Debug', or 'DebugRelease', but is set to '${CMAKE_BUILD_TYPE}'."
    )
endif()

if(
    "${CMAKE_BUILD_TYPE}" STREQUAL "Debug"
    OR "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease"
)
    set(PRISMS_PF_BUILD_DEBUG "ON")
else()
    set(PRISMS_PF_BUILD_DEBUG "OFF")
endif()

if(
    "${CMAKE_BUILD_TYPE}" STREQUAL "Release"
    OR "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease"
)
    set(PRISMS_PF_BUILD_RELEASE "ON")
else()
    set(PRISMS_PF_BUILD_RELEASE "OFF")
endif()

# Find deal.II installation
find_package(deal.II 9.6.0 QUIET HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR})
if(NOT ${deal.II_FOUND})
    message(
        FATAL_ERROR
        "\n*** Could not find a recent version of deal.II. ***\n"
        "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake "
        "or set an environment variable \"DEAL_II_DIR\" that contains a path to a "
        "recent version of deal.II."
    )
endif()

message(
    STATUS
    "Found deal.II version ${DEAL_II_PACKAGE_VERSION} at '${deal.II_DIR}'"
)

set(DEALII_INSTALL_VALID ON)

if(NOT DEAL_II_WITH_LAPACK)
    message(SEND_ERROR "\n**deal.II was built without support for lapack!\n")
    set(DEALII_INSTALL_VALID OFF)
endif()

if(NOT DEAL_II_WITH_P4EST)
    message(SEND_ERROR "\n**deal.II was built without support for p4est!\n")
    set(DEALII_INSTALL_VALID OFF)
endif()

if(NOT DEAL_II_WITH_MPI)
    message(SEND_ERROR "\n**deal.II was built without support for MPI!\n")
    set(DEALII_INSTALL_VALID OFF)
endif()

if(NOT DEALII_INSTALL_VALID)
    message(
        FATAL_ERROR
        "\nPRISMS-PF requires a deal.II installation with certain features enabled!\n"
    )
endif()

# Optional deal.II packages
message(STATUS "Using PRISMS_PF_WITH_ZLIB = '${PRISMS_PF_WITH_ZLIB}'")
if(PRISMS_PF_WITH_ZLIB)
    if(DEAL_II_WITH_ZLIB)
        message(STATUS "  Found deal.II installation with zlib")
    else()
        message(
            FATAL_ERROR
            "deal.II installation with zlib not found. Disable PRISMS_PF_WITH_ZLIB or recompile deal.II with zlib."
        )
    endif()
endif()

message(STATUS "Using PRISMS_PF_WITH_SUNDIALS = '${PRISMS_PF_WITH_SUNDIALS}'")
if(PRISMS_PF_WITH_SUNDIALS)
    if(DEAL_II_WITH_SUNDIALS)
        message(STATUS "  Found deal.II installation with SUNDIALS")
    else()
        message(
            FATAL_ERROR
            "deal.II installation with SUNDIALS not found. Disable PRISMS_PF_WITH_SUNDIALS or recompile deal.II with SUNDIALS."
        )
    endif()
endif()

# Load deal.II cached variables
deal_ii_initialize_cached_variables()

# Caliper
message(STATUS "Using PRISMS_PF_WITH_CALIPER = '${PRISMS_PF_WITH_CALIPER}'")
if(PRISMS_PF_WITH_CALIPER)
    find_package(CALIPER)
    if(${CALIPER_FOUND})
        include_directories(${CALIPER_INCLUDE_DIR})
        message(STATUS "  Caliper found at ${CALIPER_DIR}")
    else()
        message(
            FATAL_ERROR
            "Caliper not found. Disable PRISMS_PF_WITH_CALIPER or specify a hint to your installation directory with CALIPER_DIR"
        )
    endif()
endif()

# Make and ninja build options
if(CMAKE_GENERATOR MATCHES "Ninja")
    set(_make_command "$ ninja")
else()
    set(_make_command " $ make")
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
        build.ninja rules.ninja .ninja_deps .ninja_log output.txt summary.log
        *.vtk *.vtu *.pvtu
    COMMENT "distclean invoked"
)

# List of build types
if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    list(APPEND PRISMS_PF_BUILD_TYPES "Debug")
endif()

if("${CMAKE_BUILD_TYPE}" MATCHES "Release")
    list(APPEND PRISMS_PF_BUILD_TYPES "Release")
endif()
