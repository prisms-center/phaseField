#
# Find Caliper and run some checks
#

# TODO: Add some information about these variables
set(PRISMS_PF_CALIPER_VERSION "")
set(PRISMS_PF_CALIPER_DIR "")
set(PRISMS_PF_CALIPER_SOURCE "") # CALIPER_DIR or vendored

# There are two different ways to find Caliper:
#   1. Using CALIPER_DIR
#   2. As a vendored dependency
# Importantly, this is also the order in which we search.

if(NOT PRISMS_PF_FORCE_BUNDLED_CALIPER)
  find_package(
    caliper
    HINTS
      ${CALIPER_DIR}
      $ENV{CALIPER_DIR}
  )
  if(caliper_FOUND)
    prisms_pf_add_dependency_targets(DEBUG PUBLIC caliper)
    prisms_pf_add_dependency_targets(RELEASE PUBLIC caliper)
    set(PRISMS_PF_CALIPER_VERSION ${Caliper_VERSION})
    set(PRISMS_PF_CALIPER_DIR ${Caliper_DIR})
    set(PRISMS_PF_CALIPER_SOURCE "CALIPER_DIR")
    return()
  endif()
endif()

# If we can't find Caliper, vendor it
# TODO: Implement this
