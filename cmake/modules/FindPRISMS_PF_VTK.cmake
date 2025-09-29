#
# Try and find the vtk library
#

# If deal.II wasn't built with VTK we have to check for the installation
# elsewhere
if(NOT DEAL_II_WITH_VTK AND PRISMS_PF_WITH_VTK)
    message(
        STATUS
        "deal.II wasn't built with VTK, checking for separate installation"
    )

    find_package(VTK QUIET HINTS ${VTK_DIR} $ENV{VTK_DIR})
    set(VTK_VERSION "${VTK_VERSION}")
    set(VTK_MAJOR_VERSION "${VTK_MAJOR_VERSION}")
    set(VTK_MINOR_VERSION "${VTK_MINOR_VERSION}")
    set(VTK_INCLUDE_DIR
        ${VTK_PREFIX_PATH}/include/vtk-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}
    )
    # Filter the included libraries
    set(_libraries)
    foreach(_library ${VTK_LIBRARIES})
        if(
            NOT ${_library} MATCHES "Python"
            AND NOT ${_library} MATCHES "MPI4Py"
        )
            get_target_property(
                _configurations
                ${_library}
                IMPORTED_CONFIGURATIONS
            )
            if(_configurations)
                foreach(_configuration ${_configurations})
                    get_target_property(
                        _imported_location
                        ${_library}
                        IMPORTED_LOCATION_${_configuration}
                    )
                    list(APPEND _libraries ${_imported_location})
                endforeach()
            endif()
        endif()
    endforeach()
    set(VTK_NEW_LIBRARIES ${_libraries})
    set(VTK_BUILT_SEPARATELY ON)
endif()
