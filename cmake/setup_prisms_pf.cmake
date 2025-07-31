#
# Setup PRISMS-PF definitions
#

# General info
set_if_empty(PRISMS_PF_PACKAGE_NAME "PRISMS-PF")

# Grab the version
file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" _version LIMIT_COUNT 1)
set_if_empty(PRISMS_PF_VERSION "${_version}")

# Decompose the version
string(
    REGEX REPLACE
    "^([0-9]+)\\..*"
    "\\1"
    PRISMS_PF_VERSION_MAJOR
    "${PRISMS_PF_VERSION}"
)
string(
    REGEX REPLACE
    "^[0-9]+\\.([0-9]+).*"
    "\\1"
    PRISMS_PF_VERSION_MINOR
    "${PRISMS_PF_VERSION}"
)
string(
    REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.([0-9]+).*"
    "\\1"
    PRISMS_PF_VERSION_SUBMINOR
    "${PRISMS_PF_VERSION}"
)

# List of build types
if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    list(APPEND PRISMS_PF_BUILD_TYPES "Debug")
endif()

if("${CMAKE_BUILD_TYPE}" MATCHES "Release")
    list(APPEND PRISMS_PF_BUILD_TYPES "Release")
endif()

# Set the standard to C++20
set_cpp_standard(20)

# Create compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(FORCE_COLORED_OUTPUT
    ON
    CACHE BOOL
    "Forces colored output when compiling with gcc and clang."
)
