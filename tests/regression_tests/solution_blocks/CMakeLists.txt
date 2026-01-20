##
#  CMake script for the PRISMS-PF applications
##

cmake_minimum_required(VERSION 3.18)

# Create a project for the application
project(
  myapp
  C
  CXX
)

# Find the PRISMS-PF package
find_package(
  prisms_pf
  QUIET
  HINTS
    ${PRISMS_PF_DIR}
    $ENV{PRISMS_PF_DIR}
)
if(PRISMS_PF_FOUND)
  # Include setup script
  include(${PRISMS_PF_CMAKE_DIR}/setup_application.cmake)

  # Set up the application
  prisms_pf_autopilot(${PRISMS_PF_CMAKE_DIR})
else()
  message(FATAL_ERROR "PRISMS-PF package not found. Please make sure it is in path.")
endif()
