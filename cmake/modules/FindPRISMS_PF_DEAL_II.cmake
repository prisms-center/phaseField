#
# Try and find the deal.II library
#

set(DEAL_II_DIR "" CACHE PATH "An optional hint to a deal.II directory")
set_if_empty(DEAL_II_DIR "$ENV{DEAL_II_DIR}")

# Try to find deal.II
message(STATUS "Looking for deal.II installation...")
message(STATUS "  DEAL_II_DIR: ${DEAL_II_DIR}")
find_package(deal.II 9.6.0 QUIET HINTS ${DEAL_II_DIR})
if(NOT ${deal.II_FOUND})
  message(
    FATAL_ERROR
    "\n*** Could not find a recent version of deal.II. ***\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake "
    "or set an environment variable \"DEAL_II_DIR\" that contains a path to a "
    "recent version of deal.II."
  )
endif()
message(STATUS "  deal.II version ${DEAL_II_PACKAGE_VERSION} at '${deal.II_DIR}'")

# Check the deal.II build types match the ones we are requesting
string(
  FIND "${DEAL_II_BUILD_TYPE}"
  "${CMAKE_BUILD_TYPE}"
  _pos
)
if(_pos EQUAL -1)
  message(
    FATAL_ERROR
    "Mismatch between the build types that deal.II was built with and "
    "the ones being requested. PRISMS-PF is tryng to build with ${CMAKE_BUILD_TYPE}, but "
    "deal.II was built with ${DEAL_II_BUILD_TYPE}."
  )
endif()

# Required deal.II features
set(DEAL_II_INSTALL_VALID ON)

if(NOT DEAL_II_WITH_LAPACK)
  message(SEND_ERROR "\n**deal.II was built without support for lapack!\n")
  set(DEAL_II_INSTALL_VALID OFF)
endif()

if(NOT DEAL_II_WITH_P4EST)
  message(SEND_ERROR "\n**deal.II was built without support for p4est!\n")
  set(DEAL_II_INSTALL_VALID OFF)
endif()

if(NOT DEAL_II_WITH_MPI)
  message(SEND_ERROR "\n**deal.II was built without support for MPI!\n")
  set(DEAL_II_INSTALL_VALID OFF)
endif()

if(NOT DEAL_II_INSTALL_VALID)
  message(
    FATAL_ERROR
    "\nPRISMS-PF requires a deal.II installation with certain features enabled!\n"
  )
endif()

# Optional deal.II features

# zlib should be required by p4est in deal.II, so it's on by default
if(PRISMS_PF_WITH_ZLIB)
  if(DEAL_II_WITH_ZLIB)
    message(STATUS "    Found deal.II installation with zlib")
  else()
    message(
      FATAL_ERROR
      "deal.II installation with zlib not found. Disable PRISMS_PF_WITH_ZLIB or recompile deal.II with zlib."
    )
  endif()
else()
  message(
    AUTHOR_WARNING
    "PRISMS_PF_WITH_ZLIB is OFF; zlib support is supposed to be disabled, but likely isn't in the code"
  )
endif()

function(remove_std_flag variable_name)
  # Get the current value
  set(flags "${${variable_name}}")

  # Remove any -std=c++XX or -std=gnu++XX flags
  string(REGEX REPLACE "-std=[^ ]+" "" flags "${flags}")

  # Clean up any extra spaces
  string(REGEX REPLACE "  +" " " flags "${flags}")
  string(STRIP "${flags}" flags)

  # Set the modified value back to the parent scope
  set(${variable_name} "${flags}" PARENT_SCOPE)
endfunction()

# Grab relevant deal.II flags and put them into our own variables
append_flags(DEAL_II_CXX_FLAGS PRISMS_PF_CXX_FLAGS)
append_flags(DEAL_II_CXX_FLAGS_DEBUG PRISMS_PF_CXX_FLAGS_DEBUG)
append_flags(DEAL_II_CXX_FLAGS_RELEASE PRISMS_PF_CXX_FLAGS_RELEASE)

# Remove and standard flags that deal.II might have because we want to be able
# to set our standard separately.
remove_std_flag(PRISMS_PF_CXX_FLAGS)
remove_std_flag(PRISMS_PF_CXX_FLAGS_DEBUG)
remove_std_flag(PRISMS_PF_CXX_FLAGS_RELEASE)

append_flags(DEAL_II_LINKER_FLAGS PRISMS_PF_LINKER_FLAGS)
append_flags(DEAL_II_LINKER_FLAGS_DEBUG PRISMS_PF_LINKER_FLAGS_DEBUG)
append_flags(DEAL_II_LINKER_FLAGS_RELEASE PRISMS_PF_LINKER_FLAGS_RELEASE)
