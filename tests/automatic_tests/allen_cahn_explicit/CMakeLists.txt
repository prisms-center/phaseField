##
#  CMake script for the PRISMS-PF applications
#  Adapted from the ASPECT CMake file
##

cmake_minimum_required(VERSION 3.13.4)

# Create a project for the application
project(myapp CXX)

# Get the core library's directory (3 levels up)
get_filename_component(
    PRISMS_PF_CORE_DIR
    ${CMAKE_CURRENT_LIST_DIR}/../../..
    ABSOLUTE
)

# Find the PRISMS-PF package
find_package(prisms_pf REQUIRED HINTS ${PRISMS_PF_CORE_DIR})
if(PRISMS_PF_FOUND)
    # Include setup script
    include(${PRISMS_PF_CORE_DIR}/cmake/setup_application.cmake)

    # Set up the application
    prisms_pf_autopilot(${PRISMS_PF_CORE_DIR})
endif()
