##
#  CMake script for the PRISMS-PF applications
#  Adapted from the ASPECT CMake file
##

cmake_minimum_required(VERSION 3.18)

# Get the core library's directory (3 levels up from the macro file)
get_filename_component(PRISMS_PF_CORE_DIR ${CMAKE_CURRENT_LIST_DIR}/../../.. ABSOLUTE)

# Include core library's configuration
include(${PRISMS_PF_CORE_DIR}/cmake/prisms_pf_config.cmake)

# Include setup script
include(${PRISMS_PF_CORE_DIR}/cmake/setup_application.cmake)

# Create a project for the application
project(myapp CXX)

# Set the override for the src files
set(TARGET_SRC_OVERRIDE "${CMAKE_CURRENT_SOURCE_DIR}/main.cc")

# Include the autopilot macro
include(${PRISMS_PF_CORE_DIR}/cmake/macros/macro_prisms_pf_autopilot.cmake)

# Set up the application
prisms_pf_autopilot(${PRISMS_PF_CORE_DIR})
