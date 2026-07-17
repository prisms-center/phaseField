#
# Configuration options setup
#

# Decompose the version information
version_decomposition(${PRISMS_PF_VERSION})

# Grab git version information
prisms_pf_git_version()

# Some general configuration options
option(PRISMS_PF_AUTODETECTION "Autodetection of PRISMS-PF dependencies." ON)

option(PRISMS_PF_UNIT_TESTS "Whether to build the unit tests or not." OFF)
option(PRISMS_PF_REGRESSION_TESTS "Whether to build the regression tests or not." OFF)
option(PRISMS_PF_PERFORMANCE_TESTS "Whether to build the performance tests or not." OFF)
option(PRISMS_PF_EXAMPLES "Whether to build examples or not." OFF)
option(PRISMS_PF_DOCS "Whether to build documentation or not." OFF)

# CI configuration options
option(PRISMS_PF_CLANG_TIDY "Whether to run clang-tidy during the build stage." OFF)

# Backend configuration options
# NOTE: This relies on the fact that deal.II is also compiled with 64-bit indices
option(PRISMS_PF_64BIT_INDICES "Whether to use 64-bit indices for large simulations" OFF)
# NOTE: This can only go higher than 6 if deal.II is compiled with >6 for the element degree.
# TODO: Actually implement a check for deal.II having 6+ element degree
option(PRISMS_PF_ADDITIONAL_DEGREES "Whether to compile with element degrees of 3+." OFF)
# NOTE: When on deal.II will try to use the max number of threads with respect to the number
# of MPI processes.
option(PRISMS_PF_THREADS "Whether to use threading." ON)
option(PRISMS_PF_GPU "Whether to use GPU backends" OFF)
if(PRISMS_PF_GPU)
  message(WARNING "GPU backends are still under development")
endif()

# Additional dependencies
option(
  PRISMS_PF_WITH_VTK
  "Whether the user wants to compile PRISMS-PF with VTK, or not"
  OFF
)
option(PRISMS_PF_FORCE_BUNDLED_VTK "Whether the user wants to force the bundling on VTK, or not" OFF)

option(
  PRISMS_PF_WITH_HDF5
  "Whether the user wants to compile PRISMS-PF with HDF5, or not"
  OFF
)

option(
  PRISMS_PF_WITH_CALIPER
  "Whether the user wants to compile PRISMS-PF with the profiling code Caliper, or not."
  OFF
)
option(PRISMS_PF_FORCE_BUNDLED_CALIPER "Whether the user wants to force the bundling on Caliper, or not" OFF)


# With all the options out of the way we define some global variables here
# that we'll use later. These really arise from the fact that we have to
# deal with the monstrosity that is DebugRelease and how that interacts
# with our vendored packages.

# First define targets
set(PRISMS_PF_TARGETS prisms_pf)
if(CMAKE_BUILD_TYPE STREQUAL "DebugRelease")
	list(APPEND PRISMS_PF_TARGETS prisms_pf_debug)
endif()

# Next define vendored packages
set(PRISMS_PF_VENDORED_PACKAGES libassert)
if(PRISMS_PF_UNIT_TESTS OR PRISMS_PF_REGRESSION_TESTS OR PRISMS_PF_PERFORMANCE_TESTS)
  list(APPEND PRISMS_PF_VENDORED_PACKAGES Catch2)
endif()
if(PRISMS_PF_WITH_VTK AND PRISMS_PF_FORCE_BUNDLED_VTK)
	list(APPEND PRISMS_PF_VENDORED_PACKAGES VTK)
endif()
if(PRISMS_PF_WITH_CALIPER AND PRISMS_PF_FORCE_BUNDLED_CALIPER)
	list(APPEND PRISMS_PF_VENDORED_PACKAGES Caliper)
endif()

# Once we have the vendored packages, we can define packages that belong
# to the release and debug targets, if they exist. We follow the same
# naming convention as above. PRISMS_PF_PACKAGES is shared between
# the release and debug lists. The exception to this is deal.II since
# they also have a DebugRelease target. Vendored packages that follow
# the norm must be built with release and debug configurations in 
# DebugRelease
set(PRISMS_PF_PACKAGES dealii)
if(PRISMS_PF_WITH_VTK)
  list(APPEND PRISMS_PF_PACKAGES VTK)
endif()
if(PRISMS_PF_WITH_CALIPER)
	list(APPEND PRISMS_PF_PACKAGES Caliper)
endif()
set(PRISMS_PF_PACKAGES_RELEASE "")
set(PRISMS_PF_PACKAGES_DEBUG "")
foreach(pkg IN LISTS PRISMS_PF_VENDORED_PACKAGES)
	list(APPEND PRISMS_PF_PACKAGES_RELEASE ${pkg})
  list(APPEND PRISMS_PF_PACKAGES_DEBUG "${pkg}_debug")
endforeach()
list(REMOVE_ITEM PRISMS_PF_PACKAGES ${PRISMS_PF_PACKAGES_RELEASE})

# TODO: Handle these
set(PRISMS_PF_CXX_FLAGS "")
set(PRISMS_PF_CXX_FLAGS_DEBUG "")
set(PRISMS_PF_CXX_FLAGS_RELEASE "")

set(PRISMS_PF_LINKER_FLAGS "")
set(PRISMS_PF_LINKER_FLAGS_DEBUG "")
set(PRISMS_PF_LINKER_FLAGS_RELEASE "")

set(PRISMS_PF_LOG_SUMMARY "") # for summary.log

# TODO: Move this?
# Write information to log
file(REMOVE ${CMAKE_BINARY_DIR}/summary.log)
file(REMOVE ${CMAKE_BINARY_DIR}/detailed.log)

format_subsection(_tmp "prisms-pf cmake summary")
print_to_both("${_tmp}")
format_variable(_tmp PRISMS_PF_TARGETS)
print_to_both("${_tmp}")
format_variable(_tmp PRISMS_PF_VENDORED_PACKAGES)
print_to_both("${_tmp}")
format_variable(_tmp PRISMS_PF_PACKAGES)
print_to_both("${_tmp}")
format_variable(_tmp PRISMS_PF_PACKAGES_RELEASE)
print_to_both("${_tmp}")
format_variable(_tmp PRISMS_PF_PACKAGES_DEBUG)
print_to_both("${_tmp}")
