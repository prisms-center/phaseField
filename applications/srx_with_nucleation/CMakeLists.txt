##
#  CMake script for the PRISMS-PF applications
#  Adapted from the ASPECT CMake file
##

cmake_minimum_required(VERSION 3.8.0)

include(${CMAKE_SOURCE_DIR}/../../cmake/setup_application.cmake)

project(myapp CXX)

# Set location of files
include_directories(${CMAKE_SOURCE_DIR}/../../include)
include_directories(${CMAKE_SOURCE_DIR}/../../src)
include_directories(${CMAKE_SOURCE_DIR})

# Set the location of the main.cc file
set(TARGET_SRC "${CMAKE_SOURCE_DIR}/../main.cc")

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