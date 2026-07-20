#
# Find deal.II and run some checks
#

# Once we've found MPI, we need to find deal.II. Since deal.II likes to do things
# a little differently, provide some directory hints.
find_package(
  deal.II
  9.6.0
  REQUIRED
  HINTS
    ${DEAL_II_DIR}
    $ENV{DEAL_II_DIR}
)

# TODO: Don't rely on deal.II do stuff ourselves

function(strip_dealii_flags input output)
  string(REGEX REPLACE "-O[0-9s]|-std=[^ ]+" "" stripped "${input}")
  string(STRIP "${stripped}" stripped)
  set(${output} "${stripped}" PARENT_SCOPE)
endfunction()

strip_dealii_flags("${DEAL_II_CXX_FLAGS}" DEAL_II_CXX_FLAGS)
strip_dealii_flags("${DEAL_II_CXX_FLAGS_DEBUG}" DEAL_II_CXX_FLAGS_DEBUG)
strip_dealii_flags("${DEAL_II_CXX_FLAGS_RELEASE}" DEAL_II_CXX_FLAGS_RELEASE)

# Check the deal.II was built with the right dependencies
if(NOT DEAL_II_WITH_MPI)
  message(FATAL_ERROR "deal.II was not built with MPI support!")
endif()
if(NOT DEAL_II_WITH_LAPACK)
  message(FATAL_ERROR "deal.II was not built with LAPACK support!")
endif()
if(NOT DEAL_II_WITH_ZLIB)
  message(FATAL_ERROR "deal.II was not built with zlib support!")
endif()
if(NOT DEAL_II_WITH_P4EST)
  message(FATAL_ERROR "deal.II was not built with p4est support!")
endif()

# Check that deal.II links against the same MPI that we find earlier.
get_target_property(DEAL_II_MPI_LIBS dealii::interface_mpi INTERFACE_LINK_LIBRARIES)
get_target_property(PRISMS_PF_MPI_LIBS MPI::MPI_CXX INTERFACE_LINK_LIBRARIES)

if(NOT DEAL_II_MPI_LIBS STREQUAL PRISMS_PF_MPI_LIBS)
  message(
    WARNING
    "deal.II MPI libraries: ${DEAL_II_MPI_LIBS}\n"
    "Found MPI libraries:   ${PRISMS_PF_MPI_LIBS}\n"
    "These do not match! This may cause runtime issues."
  )
endif()

# Check that deal.II has the same width indices as PRISMS-PF
if(NOT DEAL_II_WITH_64BIT_INDICES STREQUAL PRISMS_PF_64BIT_INDICES)
  message(
    FATAL_ERROR
    "PRISMS-PF and deal.II must both use the same index width. You have:\n"
    "  PRISMS_PF_64BIT_INDICES = ${PRISMS_PF_64BIT_INDICES}\n"
    "  DEAL_II_WITH_64BIT_INDICES = ${DEAL_II_WITH_64BIT_INDICES}\n"
  )
endif()

# Add deal.II to the Release and Debug lists
# Again we run into some more of the uniqueness of deal.II with the way they set
# up build types. Their default build type is `DebugRelease`, which builds both
# Debug and Release targets. As such, our project much choose what deal.II target
# to link to when the user compiles with Debug and Release.
#
# When deal.II is built with only Debug or only Release, we link against for all
# PRISMS-PF builds. For example, if deal.II is built in Release, PRISMS-PF in
# Debug will only contain assertions and debugging for our library. In most cases,
# this is sufficient to catch most issues.
#
# When deal.II is build with DebugRelease, we link to the target that corresponds
# to the PRISMS-PF build type. This is recommended for developers.

if(DEAL_II_BUILD_TYPE STREQUAL "DebugRelease")
  if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    prisms_pf_add_dependency_target(dealii::dealii_debug DEBUG PUBLIC)
  elseif(${CMAKE_BUILD_TYPE} STREQUAL "DebugRelease")
    prisms_pf_add_dependency_target(dealii::dealii_debug DEBUG PUBLIC)
    prisms_pf_add_dependency_target(dealii::dealii_release RELEASE PUBLIC)
  else()
    prisms_pf_add_dependency_target(dealii::dealii_release RELEASE PUBLIC)
  endif()
else()
  prisms_pf_add_dependency_target(dealii::dealii DEBUG PUBLIC)
  prisms_pf_add_dependency_target(dealii::dealii RELEASE PUBLIC)
endif()
