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
