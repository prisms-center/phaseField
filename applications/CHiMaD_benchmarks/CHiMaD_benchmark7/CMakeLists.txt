##
#  CMake script for the PRISMS-PF applications:
##

CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)

# Find deal.II installation
FIND_PACKAGE(deal.II 9.2.0 QUIET REQUIRED
	HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

# Check to make sure deal.II is configured with p4est
if(NOT ${DEAL_II_WITH_P4EST})
  MESSAGE(FATAL_ERROR "deal.II must be installed with p4est.")
endif()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# Set up the debug, release, and run targets
ADD_CUSTOM_TARGET(debug
  COMMAND +env ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND +env ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )
ADD_CUSTOM_TARGET(release
  COMMAND +env ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND +env ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )
ADD_CUSTOM_TARGET(run COMMAND main
  COMMENT "Run with ${CMAKE_BUILD_TYPE} configuration"
  )

PROJECT(myapp)

# Add postprocess.cc and nucleation.cc if they exist
if(EXISTS "postprocess.cc")
	add_definitions(-DPOSTPROCESS_FILE_EXISTS)
endif()
if(EXISTS "nucleation.cc")
	add_definitions(-DNUCLEATION_FILE_EXISTS)
endif()

# Add main.cc for executable target
ADD_EXECUTABLE(main main.cc )

# deal.II linker
DEAL_II_SETUP_TARGET(main)

# Link libraries for the build type
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
	TARGET_LINK_LIBRARIES(main ${CMAKE_SOURCE_DIR}/../../../libprisms_pf.a)
elseif(${CMAKE_BUILD_TYPE} STREQUAL "DebugRelease")
	TARGET_LINK_LIBRARIES(main ${CMAKE_SOURCE_DIR}/../../../libprisms_pf.a)
else()
	TARGET_LINK_LIBRARIES(main ${CMAKE_SOURCE_DIR}/../../../libprisms_pf_debug.a)
endif()
