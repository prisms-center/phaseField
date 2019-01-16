##
#  CMake script for the PRISMS-PF applications:
##

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

# Find deal.II installation
FIND_PACKAGE(deal.II 8.3.0 REQUIRED
	HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

# Check to make sure deal.II is configured with p4est
IF(NOT ${DEAL_II_WITH_P4EST})
  MESSAGE(FATAL_ERROR "\n"
    "*** deal.II was not installed with p4est. ***\n\n"
    “The p4est library is a mandatory prerequisite for PRISMS-PF. Please consult the \n”
    “user guide to confirm that deal.II and p4est were installed and configured correctly.”
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

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
if (${CMAKE_BUILD_TYPE} MATCHES DebugRelease)
	SET(CMAKE_BUILD_TYPE Debug)
endif()

# Check if postprocess.cc and nucleation.cc exist and set preprocessor variables
if (EXISTS "postprocess.cc")
	add_definitions(-DPOSTPROCESS_FILE_EXISTS)
endif()
if (EXISTS "nucleation.cc")
	add_definitions(-DNUCLEATION_FILE_EXISTS)
endif()

# Append extra flags for the GNU compiler to suppress some warnings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(DEAL_II_CXX_FLAGS_DEBUG "${DEAL_II_CXX_FLAGS_DEBUG} -Wno-maybe-uninitialized -Wno-unused-parameter -Wno-extra")
        set(DEAL_II_CXX_FLAGS_RELEASE "${DEAL_II_CXX_FLAGS_RELEASE} -Wno-maybe-uninitialized -Wno-unused-parameter -Wno-extra")
        #set(DEAL_II_CXX_FLAGS_DEBUG "${DEAL_II_CXX_FLAGS_DEBUG} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")
        #set(DEAL_II_CXX_FLAGS_RELEASE "${DEAL_II_CXX_FLAGS_RELEASE} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")
endif()

ADD_EXECUTABLE(main main.cc )

DEAL_II_SETUP_TARGET(main)

set(cmd "cmake")
set(arg "CMakeLists.txt")
set(dir ${PROJECT_SOURCE_DIR}/../..)
EXECUTE_PROCESS(COMMAND ${cmd} ${arg}
		WORKING_DIRECTORY ${dir})

set(cmd "make")

EXECUTE_PROCESS(COMMAND ${cmd}
		WORKING_DIRECTORY ${dir})

if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
	TARGET_LINK_LIBRARIES(main ${CMAKE_SOURCE_DIR}/../../libprisms_pf.a)
elseif(${CMAKE_BUILD_TYPE} STREQUAL "DebugRelease")
	TARGET_LINK_LIBRARIES(main ${CMAKE_SOURCE_DIR}/../../libprisms_pf.a)
else()
	TARGET_LINK_LIBRARIES(main ${CMAKE_SOURCE_DIR}/../../libprisms_pf_debug.a)
endif()
