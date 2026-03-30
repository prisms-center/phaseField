##
#  CMake script for the PRISMS-PF applications
##

cmake_minimum_required(VERSION 3.28)
cmake_policy(VERSION 3.28)

# Create a project for the application
project(my_app CXX)

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
add_executable(main)

# Add sources
target_sources(
  main
  PRIVATE
    main.cc
    custom_pde.h
)

# Link PRISMS-PF
target_link_libraries(main PRIVATE PRISMS-PF::prisms_pf)
