##
#  CMake script for the phaseField applications:
##

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
FIND_PACKAGE(deal.II 8.3.0 REQUIRED)

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(myproject)

INCLUDE_DIRECTORIES(../../include)

FILE(GLOB headers ${CMAKE_SOURCE_DIR}/../../include/*.h)


#FILE(GLOB sources ${CMAKE_SOURCE_DIR}/../../src/matrixfree/matrixFreePDE.cc)

#ADD_LIBRARY(mylib ${CMAKE_SOURCE_DIR}/../../src/matrixfree/matrixFreePDE.cc )
ADD_LIBRARY(mylib ../../src/matrixfree/matrixFreePDE.cc ../../src/utilities/vectorLoad.cc ../../src/utilities/vectorBCFunction.cc)
#ADD_LIBRARY(mylib ../../src/matrixfree/matrixFreePDE.cc)

DEAL_II_SETUP_TARGET(mylib)



ADD_EXECUTABLE(main main.cc )

DEAL_II_SETUP_TARGET(main)

TARGET_LINK_LIBRARIES(main mylib)
