#
# Setup cached variables. These are just the configuration
# options for PRISMS-PF among some other things
#

# Basic configuration options
option(PRISMS_PF_AUTODETECTION "Autodetection of PRISMS-PF dependencies." ON)
message(STATUS "PRISMS_PF_AUTODETECTION = ${PRISMS_PF_AUTODETECTION}")

option(
    64BIT_INDICES
    "Whether to compile PRISMS-PF with 64bit numbers for large simulations"
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
    REDUCED_TEMPLATES
    "Wether the user wants to enable the reduced templates or not."
    OFF
)
message(STATUS "REDUCED_TEMPLATES = ${REDUCED_TEMPLATES}")
mark_as_advanced(REDUCED_TEMPLATES)

option(
    ADDITIONAL_DEGREES
    "Wether the user wants to enable the compilation of additional element degrees or not."
    OFF
)
if(REDUCED_TEMPLATES)
    set(ADDITIONAL_DEGREES OFF)
endif()
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

# Hide some cmake variables so we can inherit them from deal.II
set(PRISMS_PF_REMOVED_FLAGS
    CMAKE_CXX_FLAGS
    CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_MINSIZEREL
    CMAKE_CXX_FLAGS_RELWITHDEBINFO
    CMAKE_SHARED_LINKER_FLAGS
    CMAKE_SHARED_LINKER_FLAGS_DEBUG
    CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL
    CMAKE_SHARED_LINKER_FLAGS_RELEASE
    CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO
)
foreach(_flag ${PRISMS_PF_REMOVED_FLAGS})
    set(${_flag} ${${_flag}} CACHE INTERNAL "" FORCE)
endforeach()

# Cache some cmake variables of ours and mark them as advanced
set(PRISMS_PF_FLAGS
    PRISMS_PF_CXX_FLAGS
    PRISMS_PF_CXX_FLAGS_DEBUG
    PRISMS_PF_CXX_FLAGS_RELEASE
    PRISMS_PF_LINKER_FLAGS
    PRISMS_PF_LINKER_FLAGS_DEBUG
    PRISMS_PF_LINKER_FLAGS_RELEASE
)
foreach(_flag ${PRISMS_PF_FLAGS})
    set(${_flag}
        "${${_flag}}"
        CACHE STRING
        "These user-supplied variables are appended to those inherited by deal.II"
    )
    mark_as_advanced(${_flag})
endforeach()

# Set compiler flags into the respective PRISMS-PF variable
foreach(_flag CXX_FLAG CXX_FLAGS_DEBUG CXX_FLAGS_RELEASE)
    if(NOT "${CMAKE_SHARED_${_flag}}" STREQUAL "")
        message(
            STATUS
            "Prepending \${CMAKE_SHARED_${_flag}} to \${PRISMS_PF_${_flag}}"
        )
        set(PRISMS_PF_${_flag} "${CMAKE_${_flag}} ${PRISMS_PF_${_flag}}")
    endif()
endforeach()

# Do the same for linker flags
foreach(_flag LINKER_FLAGS LINKER_FLAGS_DEBUG LINKER_FLAGS_RELEASE)
    if(NOT "${CMAKE_SHARED_${_flag}}" STREQUAL "")
        message(
            STATUS
            "Prepending \${CMAKE_SHARED_${_flag}} to \${PRISMS_PF_${_flag}}"
        )
        set(PRISMS_PF_${_flag} "${CMAKE_${_flag}} ${PRISMS_PF_${_flag}}")
    endif()
endforeach()

# Store these flags into a SAVED variable so we can
# cross-reference with the deal.II ones that we inherit
# later. Also, set to empty string.
foreach(_flag ${PRISMS_PF_FLAGS})
    set(${_flag}_SAVED ${${_flag}})
    set(${_flag} "")
endforeach()

# Also set the removed flags to empty
foreach(_flag ${PRISMS_PF_REMOVED_FLAGS})
    set(${_flag} "")
endforeach()

# Read the environmental flags and add them to the saved variables
set(PRISMS_PF_CXX_FLAGS_SAVED "$ENV{CXXFLAGS} ${PRISMS_PF_CXX_FLAGS_SAVED}")
string(STRIP "${PRISMS_PF_CXX_FLAGS_SAVED}" PRISMS_PF_CXX_FLAGS_SAVED)
set(PRISMS_PF_LINKER_FLAGS_SAVED
    "$ENV{LDFLAGS} ${PRISMS_PF_LINKER_FLAGS_SAVED}"
)
string(STRIP "${PRISMS_PF_LINKER_FLAGS_SAVED}" PRISMS_PF_LINKER_FLAGS_SAVED)
unset(ENV{CXXFLAGS})
unset(ENV{LDFLAGS})
unset(ENV{NVCCFLAGS})
