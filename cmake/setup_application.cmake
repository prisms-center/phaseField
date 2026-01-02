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

# Make and ninja build options
if(CMAKE_GENERATOR MATCHES "Ninja")
  set(_make_command "$ ninja")
else()
  set(_make_command " $ make")
endif()

# Setup the build type (debug, release, debugrelease)
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(
    CMAKE_BUILD_TYPE
    "DebugRelease"
    CACHE STRING
    "Choose the type of build, options are: Debug, Release and DebugRelease."
    FORCE
  )
endif()

# List of build types for the application
if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  list(APPEND APPLICATION_BUILD_TYPES "Debug")
endif()

if("${CMAKE_BUILD_TYPE}" MATCHES "Release")
  list(APPEND APPLICATION_BUILD_TYPES "Release")
endif()

# Check that the core library was compiled with those build types
foreach(BUILD_TYPE ${APPLICATION_BUILD_TYPES})
  if(NOT BUILD_TYPE IN_LIST PRISMS_PF_BUILD_TYPES)
    message(FATAL_ERROR "Build type '${BUILD_TYPE}' not found in PRISMS_PF_BUILD_TYPES")
  endif()
endforeach()

# For each of the build types as some targets that switch the type and compile
if("Debug" IN_LIST APPLICATION_BUILD_TYPES)
  add_custom_target(
    debug
    COMMAND
      ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND
      ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Configuring and compiling in DEBUG mode..."
  )
  if("Release" IN_LIST APPLICATION_BUILD_TYPES)
    add_custom_target(
      debugrelease
      COMMAND
        ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
      COMMAND
        ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMENT "Configuring and compiling in DEBUG + RELEASE mode..."
    )
  endif()
endif()
if("Release" IN_LIST APPLICATION_BUILD_TYPES)
  add_custom_target(
    release
    COMMAND
      ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=DebugRelease ${CMAKE_SOURCE_DIR}
    COMMAND
      ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Configuring and compiling in RELEASE mode..."
  )
endif()
