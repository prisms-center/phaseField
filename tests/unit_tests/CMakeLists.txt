##
#  CMake for the PRISMS-PF unit tests
#  Adapted from the ASPECT CMake file
##

cmake_minimum_required(VERSION 3.1.0)

project(prisms_pf_unit_tests)

# =========================================================
# Some basic bookkeeping
# =========================================================

# Check that a prior CMakeCache is not located in the build directory
if(EXISTS ${CMAKE_SOURCE_DIR}/CMakeCache.txt)
  message(FATAL_ERROR "Detected the file:\n"
  "${CMAKE_SOURCE_DIR}/CMakeCache.txt\n"
  "in your source directory, which may be leftover from prior builds. " 
  "Please delete the file before running cmake again.")
endif()

# =========================================================
# External libraries
# =========================================================

# Find deal.II installation
find_package(deal.II 9.2.0 QUIET
	HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR}
  )
if(NOT ${deal.II_FOUND})
  message(FATAL_ERROR "\n*** Could not find a recent version of deal.II. ***\n"
  "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake "
  "or set an environment variable \"DEAL_II_DIR\" that contains a path to a "
  "recent version of deal.II."
  )
endif()

message(STATUS "Found deal.II version ${DEAL_II_PACKAGE_VERSION} at '${deal.II_DIR}'")

set(DEALII_INSTALL_VALID ON)

if(NOT DEAL_II_WITH_P4EST)
    message(SEND_ERROR
      "\n**deal.II was built without support for p4est!\n"
      )
    set(DEALII_INSTALL_VALID OFF)
endif()

if(NOT DEALII_INSTALL_VALID)
  message(FATAL_ERROR
    "\nPRISMS-PF requires a deal.II installation with certain features enabled!\n"
    )
endif()

deal_ii_initialize_cached_variables()

# Create compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(FORCE_COLORED_OUTPUT ON CACHE BOOL "Forces colored output when compiling with gcc and clang.")

# deal.II versions >=9.5 disable deprecation warnings in user code. Enable
# the warnings again by removing the flag that disables them.
string(REPLACE "-Wno-deprecated-declarations" "" DEAL_II_WARNING_FLAGS "${DEAL_II_WARNING_FLAGS}")

# Set additional compiler flags
set(PRISMS_PF_ADDITIONAL_CXX_FLAGS "" CACHE STRING "Additional CMAKE_CXX_FLAGS applied after the deal.II options.")

if(NOT PRISMS_PF_ADDITIONAL_CXX_FLAGS STREQUAL "")
  message(STATUS "Appending PRISMS_PF_ADDITIONAL_CXX_FLAGS: '${PRISMS_PF_ADDITIONAL_CXX_FLAGS}':")
  string(APPEND DEAL_II_CXX_FLAGS_DEBUG " ${PRISMS_PF_ADDITIONAL_CXX_FLAGS}")
  string(APPEND DEAL_II_CXX_FLAGS_RELEASE " ${PRISMS_PF_ADDITIONAL_CXX_FLAGS}")
  message(STATUS "  DEAL_II_WARNING_FLAGS: ${DEAL_II_WARNING_FLAGS}")
  message(STATUS "  DEAL_II_CXX_FLAGS_DEBUG: ${DEAL_II_CXX_FLAGS_DEBUG}")
  message(STATUS "  DEAL_II_CXX_FLAGS_RELEASE: ${DEAL_II_CXX_FLAGS_RELEASE}")
endif()

# Set the name of the project and target:
set(TARGET "main")

# Declare all source files the target consists of:
set(TARGET_SRC
  test.cc
  # You can specify additional files here!
  )

deal_II_invoke_autopilot()
