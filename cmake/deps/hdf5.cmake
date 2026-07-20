#
# Find HDF5 and run some checks
#

# Currently, we only support HDF5 through deal.II
if(NOT DEAL_II_WITH_HDF5 AND PRISMS_PF_WITH_HDF5)
  message(
    FATAL_ERROR
    "PRISMS_PF_WITH_HDF5 is enabled, but deal.II was not built with HDF5 support. "
    "Rebuild deal.II with HDF5 enabled or disable PRISMS_PF_WITH_HDF5."
  )
endif()
