##
#  CMake script for the PRISMS-PF applications
##

cmake_minimum_required(VERSION 3.28)
cmake_policy(VERSION 3.28)

# Name the application
set(APPLICATION_NAME "allen_cahn_conserved")

# Create a project for the application
project(${APPLICATION_NAME} CXX)

# Set the build type to debug is not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type" FORCE)
endif()

# Export compile commands for IDEs
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find the PRISMS-PF package
find_package(
  PRISMS-PF
  REQUIRED
  HINTS
    ${PRISMS_PF_DIR}
    $ENV{PRISMS_PF_DIR}
)

# Declare the application executable
add_executable(${APPLICATION_NAME})

# Rename the output executable to main
set_target_properties(
  ${APPLICATION_NAME}
  PROPERTIES
    OUTPUT_NAME
      main
)

# Add sources
target_sources(
  ${APPLICATION_NAME}
  PRIVATE
    main.cc
    custom_pde.h
)

# Require C++20
# TODO: Remove this
target_compile_features(${APPLICATION_NAME} PRIVATE cxx_std_20)
set_target_properties(
  ${APPLICATION_NAME}
  PROPERTIES
    CXX_STANDARD
      20
    CXX_STANDARD_REQUIRED
      ON
    CXX_EXTENSIONS
      OFF
)

# Link PRISMS-PF
# If we're in Debug configuration and PRISMS_PF_BUILD_TYPE is DebugRelease, link to the prisms_pf_debug target,
# which is only available in prisms_pf DebugRelease configuration.
# Otherwise link to the prisms_pf target, which may be debug or release
target_link_libraries(
  ${APPLICATION_NAME}
  PRIVATE
    $<IF:$<AND:$<CONFIG:Debug>,$<STREQUAL:${PRISMS_PF_BUILD_TYPE},DebugRelease>>,PRISMS-PF::prisms_pf_debug,
    PRISMS-PF::prisms_pf>
)
