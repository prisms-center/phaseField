##
#  CMake script for the PRISMS-PF applications:
##

cmake_minimum_required(VERSION 3.1.0)
project(myapp)

# Find deal.II installation
find_package(deal.II 9.2.0 QUIET REQUIRED
	HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

# Check to make sure deal.II is configured with p4est
if(NOT ${DEAL_II_WITH_P4EST})
  message(FATAL_ERROR "deal.II must be installed with p4est.")
endif()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# Set up the debug, release, and run targets
add_custom_target(debug
  COMMAND +env ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND +env ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )
add_custom_target(release
  COMMAND +env ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND +env ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )
add_custom_target(run COMMAND main
  COMMENT "Run with ${CMAKE_BUILD_TYPE} configuration"
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
set(MAIN "${CMAKE_SOURCE_DIR}/../main.cc")

# Add main.cc for executable target
add_executable(main ${MAIN})

# deal.II linker
DEAL_II_SETUP_TARGET(main)

# Link libraries for the build type
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
	target_link_libraries(main ${CMAKE_SOURCE_DIR}/../../libprisms_pf.a)
elseif(${CMAKE_BUILD_TYPE} STREQUAL "DebugRelease")
  target_link_libraries(main ${CMAKE_SOURCE_DIR}/../../libprisms_pf.a)
else()
  target_link_libraries(main ${CMAKE_SOURCE_DIR}/../../libprisms_pf_debug.a)
endif()
