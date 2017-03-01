##
#  CMake script for the phaseField applications:
##

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
FIND_PACKAGE(deal.II 8.3.0 REQUIRED)

DEAL_II_INITIALIZE_CACHED_VARIABLES()

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )
  
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(DEAL_II_CXX_FLAGS_DEBUG "${DEAL_II_CXX_FLAGS_DEBUG} -Wno-maybe-uninitialized -Wno-unused-parameter")
        set(DEAL_II_CXX_FLAGS_RELEASE "${DEAL_II_CXX_FLAGS_RELEASE} -Wno-maybe-uninitialized -Wno-unused-parameter")
        #set(DEAL_II_CXX_FLAGS_DEBUG "${DEAL_II_CXX_FLAGS_DEBUG} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")
        #set(DEAL_II_CXX_FLAGS_RELEASE "${DEAL_II_CXX_FLAGS_RELEASE} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")
endif()

  
PROJECT(myproject)

INCLUDE_DIRECTORIES(../../include)

FILE(GLOB headers ${CMAKE_SOURCE_DIR}/../../include/*.h)

FILE(GLOB sources ${CMAKE_SOURCE_DIR}/../../src/matrixfree/*.cc)

#ADD_LIBRARY(mylib ${CMAKE_SOURCE_DIR}/../../src/matrixfree/matrixFreePDE.cc )
#ADD_LIBRARY(mylib ../../src/matrixfree/matrixFreePDE.cc ../../src/utilities/vectorLoad.cc ../../src/utilities/vectorBCFunction.cc)
#ADD_LIBRARY(mylib ../../src/matrixfree/matrixFreePDE.cc)

ADD_LIBRARY(mylib ${sources} ../../src/utilities/vectorLoad.cc ../../src/utilities/vectorBCFunction.cc)

DEAL_II_SETUP_TARGET(mylib)

ADD_EXECUTABLE(main main.cc )

DEAL_II_SETUP_TARGET(main)

TARGET_LINK_LIBRARIES(main mylib)
