#
# Find VTK and run some checks
#

# TODO: Add some information about these variables
set(PRISMS_PF_VTK_VERSION "")
set(PRISMS_PF_VTK_DIR "")
set(PRISMS_PF_VTK_SOURCE "") # deal.II, VTK_DIR, or vendored

# There are three different ways to find VTK:
#   1. Through deal.II
#   2. Using VTK_DIR
#   3. As a vendored dependency
# Importantly, this is also the order in which we search.

# First, see if deal.II has VTK. If it does, there's nothing else we have
# to do because it already handles the propagation of dependencies.
# If the user has tried to force the bundling of VTK through an error.
if(DEAL_II_WITH_VTK)
  if(PRISMS_PF_FORCE_BUNDLED_VTK)
    message(
      FATAL_ERROR
      "Attempting to use bundled VTK; however, deal.II was built with VTK!\n"
      "Either set PRISMS_PF_FORCE_BUNDLED_VTK=OFF or recompile deal.II with"
      " DEAL_II_WITH_VTK=OFF."
    )
  endif()
  set(PRISMS_PF_VTK_VERSION ${DEAL_II_VTK_VERSION})
  set(PRISMS_PF_VTK_DIR ${VTK_DIR}) # TODO: Not sure if this is working
  set(PRISMS_PF_VTK_SOURCE "deal.II")
  return()
endif()

# Second, try and find VTK normally
if(NOT PRISMS_PF_FORCE_BUNDLED_VTK)
  find_package(
    VTK
    HINTS
      ${VTK_DIR}
      $ENV{VTK_DIR}
    COMPONENTS
      CommonCore
      IOGeometry
      IOLegacy
      IOXML
  )
  if(VTK_FOUND)
    prisms_pf_add_dependency_targets(
      DEBUG
      PUBLIC
      VTK::CommonCore
      VTK::IOGeometry
      VTK::IOLegacy
      VTK::IOXML
    )
    prisms_pf_add_dependency_targets(
      RELEASE
      PUBLIC
      VTK::CommonCore
      VTK::IOGeometry
      VTK::IOLegacy
      VTK::IOXML
    )
    set(PRISMS_PF_VTK_VERSION ${VTK_VERSION})
    set(PRISMS_PF_VTK_DIR ${VTK_DIR})
    set(PRISMS_PF_VTK_SOURCE "VTK_DIR")
    return()
  endif()
endif()

# If we can't find VTK, vendor it
# TODO: Implement this
