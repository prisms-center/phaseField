##
#  CMake script for the PRISMS-PF applications
#  Adapted from the ASPECT CMake file
##

cmake_minimum_required(VERSION 3.1.0)

project(myapp)

message(STATUS "")
message(STATUS "=========================================================")
message(STATUS "Configuring PRISMS-PF application")
message(STATUS "=========================================================")
message(STATUS "")

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
  message(STATUS "Setting up PRISMS-PF application for ${CMAKE_BUILD_TYPE} mode.")
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

# Add postprocess.cc and nucleation.cc if they exist
if(EXISTS "postprocess.cc")
	add_definitions(-DPOSTPROCESS_FILE_EXISTS)
endif()
if(EXISTS "nucleation.cc")
	add_definitions(-DNUCLEATION_FILE_EXISTS)
endif()

# Set location of files
include_directories(${CMAKE_SOURCE_DIR}/../../include)
include_directories(${CMAKE_SOURCE_DIR}/../../src)
include_directories(${CMAKE_SOURCE_DIR})

# Set the location of the main.cc file
set(TARGET_SRC "${CMAKE_SOURCE_DIR}/../main.cc")

# Check if there has been updates to main library
set(dir ${PROJECT_SOURCE_DIR}/../..)
EXECUTE_PROCESS(COMMAND "rm" "CMakeCache.txt" WORKING_DIRECTORY ${dir})
EXECUTE_PROCESS(COMMAND "cmake" "CMakeLists.txt" WORKING_DIRECTORY ${dir})
EXECUTE_PROCESS(COMMAND "make" WORKING_DIRECTORY ${dir})

# Set targets & link libraries for the build type
if(${PRISMS_PF_BUILD_DEBUG} STREQUAL "ON")
  add_executable(main_debug ${TARGET_SRC})
  set_property(TARGET main_debug PROPERTY OUTPUT_NAME main-debug)
  deal_ii_setup_target(main_debug DEBUG)
  target_link_libraries(main_debug ${CMAKE_SOURCE_DIR}/../../libprisms-pf-debug.a)
endif()

if(${PRISMS_PF_BUILD_RELEASE} STREQUAL "ON")
  add_executable(main_release ${TARGET_SRC})
  set_property(TARGET main_release PROPERTY OUTPUT_NAME main)
  deal_ii_setup_target(main_release RELEASE)
  target_link_libraries(main_release ${CMAKE_SOURCE_DIR}/../../libprisms-pf-release.a)
endif()