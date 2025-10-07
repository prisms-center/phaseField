#
# Setup cached variables. These are just the configuration
# options for PRISMS-PF among some other things
#

# Basic configuration options
option(PRISMS_PF_AUTODETECTION "Autodetection of PRISMS-PF dependencies." ON)
message(STATUS "PRISMS_PF_AUTODETECTION = ${PRISMS_PF_AUTODETECTION}")

option(UNIT_TESTS "Whether to build the unit tests or not." OFF)
message(STATUS "UNIT_TESTS = ${UNIT_TESTS}")

option(REGRESSION_TESTS "Whether to build the regression tests or not." OFF)
message(STATUS "REGRESSION_TESTS = ${REGRESSION_TESTS}")

option(
    64BIT_INDICES
    "Whether to compile PRISMS-PF with 64-bit numbers for large simulations"
    OFF
)
message(STATUS "64BIT_INDICES = ${64BIT_INDICES}")

option(
    ADDITIONAL_OPTIMIZATIONS
    "Whether the user wants to enable additional optimizations, or not."
    OFF
)
message(STATUS "ADDITIONAL_OPTIMIZATIONS = ${ADDITIONAL_OPTIMIZATIONS}")

option(
    ADDITIONAL_DEGREES
    "Wether the user wants to enable the compilation of additional element degrees or not."
    OFF
)
message(STATUS "ADDITIONAL_DEGREES = ${ADDITIONAL_DEGREES}")

option(
    UNWRAP_COMPILER
    "Whether the user wants to unwrap the compiler in compile_commands.json, or not."
    OFF
)
mark_as_advanced(UNWRAP_COMPILER)

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
    message(STATUS "Setting up PRISMS-PF for ${CMAKE_BUILD_TYPE} mode.")
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
