##
#  CMake script for the PRISMS-PF applications:
##

CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)

# Find all application directories
file(GLOB APP_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/*)

# Add subdirectories for each application, noting mixed pathes
foreach(APP_DIR ${APP_DIRS})
  if(IS_DIRECTORY ${APP_DIR})
    # Search for customPDE.h
    file(GLOB_RECURSE CUSTOM_PDE ${APP_DIR}/customPDE.h)

    # Add the directory
    add_subdirectory(${APP_DIR})

  endif(IS_DIRECTORY ${APP_DIR})
endforeach(APP_DIR ${APP_DIRS})