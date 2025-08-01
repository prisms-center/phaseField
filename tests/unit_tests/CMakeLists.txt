##
#  CMake for the PRISMS-PF unit tests
#  Adapted from the ASPECT CMake file
##

cmake_minimum_required(VERSION 3.8.0)

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

# Set the standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# =========================================================
# External libraries
# =========================================================

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

if(NOT DEAL_II_WITH_P4EST)
    message(SEND_ERROR "\n**deal.II was built without support for p4est!\n")
    set(DEALII_INSTALL_VALID OFF)
endif()

if(NOT DEALII_INSTALL_VALID)
    message(
        FATAL_ERROR
        "\nPRISMS-PF requires a deal.II installation with certain features enabled!\n"
    )
endif()

deal_ii_initialize_cached_variables()

project(prisms_pf_unit_tests)

# Check that deal.II was built with vtk or we can find the package ourselves
set(VTK_BUILT_SEPARATELY
    OFF
    CACHE BOOL
    "Whether the installed VTK library was built outside of deal.II."
)
if(NOT DEAL_II_WITH_VTK)
    find_package(VTK QUIET HINTS ${VTK_DIR} $ENV{VTK_DIR})
    if(NOT VTK_FOUND)
        message(SEND_ERROR "\n**deal.II was built without support for VTK!\n")
        set(DEALII_INSTALL_VALID OFF)
    endif()
    set(VTK_VERSION "${VTK_VERSION}")
    set(VTK_MAJOR_VERSION "${VTK_MAJOR_VERSION}")
    set(VTK_MINOR_VERSION "${VTK_MINOR_VERSION}")
    set(VTK_INCLUDE_DIR
        ${VTK_PREFIX_PATH}/include/vtk-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}
    )
    # Filter the included libraries
    set(_libraries)
    foreach(_library ${VTK_LIBRARIES})
        if(
            NOT ${_library} MATCHES "Python"
            AND NOT ${_library} MATCHES "MPI4Py"
        )
            get_target_property(
                _configurations
                ${_library}
                IMPORTED_CONFIGURATIONS
            )
            if(_configurations)
                foreach(_configuration ${_configurations})
                    get_target_property(
                        _imported_location
                        ${_library}
                        IMPORTED_LOCATION_${_configuration}
                    )
                    list(APPEND _libraries ${_imported_location})
                endforeach()
            endif()
        endif()
    endforeach()
    set(VTK_NEW_LIBRARIES ${_libraries})
    set(VTK_BUILT_SEPARATELY ON)
endif()

# Create compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(FORCE_COLORED_OUTPUT
    ON
    CACHE BOOL
    "Forces colored output when compiling with gcc and clang."
)

# Setting up tests
message(STATUS "Setting up tests with CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}")

enable_testing()

add_custom_target(tests)

# diff tool
find_program(DIFF_EXECUTABLE NAMES diff HINTS ${DIFF_DIR} PATH_SUFFIXES bin)
find_program(
    NUMDIFF_EXECUTABLE
    NAMES numdiff
    HINTS ${NUMDIFF_DIR}
    PATH_SUFFIXES bin
)

mark_as_advanced(DIFF_EXECUTABLE NUMDIFF_EXECUTABLE)

if("${TEST_DIFF}" STREQUAL "")
    if(NOT NUMDIFF_EXECUTABLE MATCHES "-NOTFOUND")
        set(TEST_DIFF ${NUMDIFF_EXECUTABLE})
    elseif(NOT DIFF_EXECUTABLE MATCHES "-NOTFOUND")
        set(TEST_DIFF ${DIFF_EXECUTABLE})
    else()
        message(
            FATAL_ERROR
            "Could not find diff or numdiff. One of those must be installed for running the testsuite./n"
            "Please specify TEST_DIFF by hand."
        )
    endif()
endif()

# Set the name of the project and target:
set(TARGET "main")

file(GLOB_RECURSE TEST_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.cc")

# Set location of Catch
include_directories(${CMAKE_SOURCE_DIR}/../../contrib/catch/)

# Set location of core library files
include_directories(${CMAKE_SOURCE_DIR}/../../include)

# Declare all source files the target consists of:
set(TARGET_SRC main.cc ${TEST_SOURCES})

# Set targets & link libraries for the build type
add_executable(main ${TARGET_SRC})
set_property(TARGET main PROPERTY OUTPUT_NAME main)
deal_ii_setup_target(main DEBUG)
target_link_libraries(main ${CMAKE_SOURCE_DIR}/../../libprisms-pf-debug.a)
if(${VTK_BUILT_SEPARATELY})
    include_directories(SYSTEM ${VTK_INCLUDE_DIR})
    target_link_libraries(main ${VTK_NEW_LIBRARIES})
endif()

# CTest
add_test(NAME PRISMS_PF_Testsuite COMMAND main)
