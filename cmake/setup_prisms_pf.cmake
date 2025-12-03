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

# Features
option(
    PRISMS_PF_WITH_ZLIB
    "Whether the user wants to compile PRISMS-PF with deal.II's zlib dependency, or not."
    ON
)

option(
    PRISMS_PF_WITH_VTK
    "Whether the user wants to compiler PRISMS-PF with vtk, or not"
    OFF
)

option(
    PRISMS_PF_WITH_CALIPER
    "Whether the user wants to compile PRISMS-PF with the profiling code Caliper, or not."
    OFF
)

# Set the compiler and linker flags
set(PRISMS_PF_WARNING_FLAGS "")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU") # GCC
    set(PRISMS_PF_WARNING_FLAGS
        -Wall
        -Wextra
        -Wpedantic
				# -Wconversion # TODO: Renable this
        # -Wsign-conversion # This is disabled because deal.II uses int when it should be uint so it just produces a lot of noise for limited utility 
        -Wshadow
        -Wnon-virtual-dtor
        -Wold-style-cast
        -Wcast-align
        -Wunused
        -Woverloaded-virtual
        -Wnull-dereference
        -Wdouble-promotion
        -Wformat=2
        -Wimplicit-fallthrough
        -Wmisleading-indentation
        -Wduplicated-cond
        -Wduplicated-branches
        -Wlogical-op
        -Wuseless-cast
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang") # Clang/AppleClang
    set(PRISMS_PF_WARNING_FLAGS
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion
        # -Wsign-conversion # This is disabled because deal.II uses int when it should be uint so it just produces a lot of noise for limited utility 
        -Wshadow
        -Wnon-virtual-dtor
        -Wold-style-cast
        -Wcast-align
        -Wunused
        -Woverloaded-virtual
        -Wnull-dereference
        -Wdouble-promotion
        -Wformat=2
        -Wimplicit-fallthrough
        -Wdocumentation
        -Winconsistent-missing-destructor-override-attribute
        -Wunreachable-code
        -Wmove
        -Wloop-analysis
        -Wcomma
        -Wrange-loop-analysis
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel") # Intel classic (icc/icpc)
    set(PRISMS_PF_WARNING_FLAGS
        -Wall
        -Wextra
        -Wcheck
        -Wshadow
        -Wunused-variable
        -Wuninitialized
        -Wremarks
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM") # Newer Intel oneAPI compilers (icx/icpx)
    set(PRISMS_PF_WARNING_FLAGS
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion
        # -Wsign-conversion # This is disabled because deal.II uses int when it should be uint so it just produces a lot of noise for limited utility 
        -Wshadow
        -Wnon-virtual-dtor
        -Wold-style-cast
        -Wcast-align
        -Wunused
        -Woverloaded-virtual
        -Wnull-dereference
        -Wdouble-promotion
        -Wformat=2
        -Wimplicit-fallthrough
    )
endif()

set(PRISMS_PF_CXX_FLAGS "")
set(PRISMS_PF_CXX_FLAGS_DEBUG "")
set(PRISMS_PF_CXX_FLAGS_RELEASE "")

set(PRISMS_PF_LINKER_FLAGS "")
set(PRISMS_PF_LINKER_FLAGS_DEBUG "")
set(PRISMS_PF_LINKER_FLAGS_RELEASE "")
