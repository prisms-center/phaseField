#
# Try and find the vtk library
#

# If deal.II wasn't built with VTK we have to check for the installation
# elsewhere
if(NOT DEAL_II_WITH_VTK)
    message(
        STATUS
        "deal.II wasn't built with VTK, checking for separate installation"
    )
endif()
