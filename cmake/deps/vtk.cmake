#
# Find VTK and run some checks
#

# Currently, we only support VTK through deal.II
if(NOT DEAL_II_WITH_VTK AND PRISMS_PF_WITH_VTK)
  message(
    FATAL_ERROR
    "PRISMS_PF_WITH_VTK is enabled, but deal.II was not built with VTK support. "
    "Rebuild deal.II with VTK enabled or disable PRISMS_PF_WITH_VTK."
  )
endif()
