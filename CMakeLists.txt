##
#  CMake script for the PRISMS-PF applications:
##

CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)

# Find deal.II installation
find_package(deal.II 9.2.0 QUIET REQUIRED
	HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})
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


# Set up the paths to the library files
include_directories(include)

file(GLOB headers include/*.h)
file(GLOB matrixfree_sources src/matrixfree/*.cc)
file(GLOB userinputparameters_sources src/userInputParameters/*.cc)
file(GLOB solverparameters_sources src/SolverParameters/*.cc)
file(GLOB equationdependencyparser_sources src/EquationDependencyParser/*.cc)
file(GLOB floodfiller_sources src/FloodFiller/*.cc)
file(GLOB orderparameterremapper_sources src/OrderParameterRemapper/*.cc)
file(GLOB simplifiedgrainrepresentation_sources src/SimplifiedGrainRepresentation/*.cc)

# Add main project
PROJECT(prisms_pf)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Append extra flags for the GNU compiler to suppress some warnings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(DEAL_II_CXX_FLAGS_DEBUG "${DEAL_II_CXX_FLAGS_DEBUG} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")
  set(DEAL_II_CXX_FLAGS_RELEASE "${DEAL_II_CXX_FLAGS_RELEASE} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")
endif()

ADD_LIBRARY(${PROJECT_NAME} ${matrixfree_sources} ${userinputparameters_sources} ${solverparameters_sources} ${equationdependencyparser_sources} ${floodfiller_sources} ${orderparameterremapper_sources} ${simplifiedgrainrepresentation_sources} src/utilities/sortIndexEntryPairList.cc src/variableAttributeLoader/variableAttributeLoader.cc src/utilities/vectorBCFunction.cc src/inputFileReader/inputFileReader.cc src/parallelNucleationList/parallelNucleationList.cc src/variableContainer/variableContainer.cc)
DEAL_II_SETUP_TARGET(${PROJECT_NAME})

PROJECT(prisms_pf_debug)
SET(CMAKE_BUILD_TYPE Debug)

# Append extra flags for the GNU compiler to suppress some warnings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(DEAL_II_CXX_FLAGS_DEBUG "${DEAL_II_CXX_FLAGS_DEBUG} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")
  set(DEAL_II_CXX_FLAGS_RELEASE "${DEAL_II_CXX_FLAGS_RELEASE} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")
endif()

ADD_LIBRARY(${PROJECT_NAME} ${matrixfree_sources} ${userinputparameters_sources} ${solverparameters_sources} ${equationdependencyparser_sources} ${floodfiller_sources} ${orderparameterremapper_sources} ${simplifiedgrainrepresentation_sources} src/utilities/sortIndexEntryPairList.cc src/variableAttributeLoader/variableAttributeLoader.cc src/utilities/vectorBCFunction.cc src/inputFileReader/inputFileReader.cc src/parallelNucleationList/parallelNucleationList.cc src/variableContainer/variableContainer.cc)
DEAL_II_SETUP_TARGET(${PROJECT_NAME})
