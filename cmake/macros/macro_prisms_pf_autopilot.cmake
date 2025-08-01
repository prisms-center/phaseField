#
# Autopilot function that will set up an application
#

macro(prisms_pf_autopilot PRISMS_PF_CORE_DIR)
    # Add the script files
    add_subdirectory(
        "${PRISMS_PF_CORE_DIR}/cmake/scripts"
        "${CMAKE_BINARY_DIR}/scripts_build"
    )

    # Enable compile commands export
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

    include_directories(${CMAKE_BINARY_DIR})

    # Set location of files
    include_directories(${PRISMS_PF_CORE_DIR}/include)
    include_directories(${PRISMS_PF_CORE_DIR}/src)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})

    # Find Caliper if enabled
    if(${PRISMS_PF_WITH_CALIPER})
        find_package(CALIPER REQUIRED)
        include_directories(${CALIPER_INCLUDE_DIR})
    endif()

    # Check that deal.II was built with vtk or we can find the package ourselves
    set(VTK_BUILT_SEPARATELY
        OFF
        CACHE BOOL
        "Whether the installed VTK library was built outside of deal.II."
    )
    if(NOT DEAL_II_WITH_VTK)
        find_package(VTK QUIET HINTS ${VTK_DIR} $ENV{VTK_DIR})
        if(NOT VTK_FOUND)
            message(
                SEND_ERROR
                "\n**deal.II was built without support for VTK!\n"
            )
            set(DEALII_INSTALL_VALID OFF)
        endif()
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

    # Set the location of application source files
    # Check if an override is specified otherwise
    # proceed with defaults
    if(NOT DEFINED TARGET_SRC_OVERRIDE)
        set(TARGET_SRC
            "${CMAKE_CURRENT_SOURCE_DIR}/main.cc"
            "${CMAKE_CURRENT_SOURCE_DIR}/equations.cc"
            "${CMAKE_CURRENT_SOURCE_DIR}/ICs_and_BCs.cc"
        )
    else()
        set(TARGET_SRC ${TARGET_SRC_OVERRIDE})
    endif()

    foreach(_build ${PRISMS_PF_BUILD_TYPES})
        string(TOLOWER ${_build} _build_lowercase)
        add_executable(main_${_build_lowercase} ${TARGET_SRC})
    endforeach()

    expand_template_instantiations(main "custom_pde.inst.in")

    # Set targets & link libraries for the build type
    if(${PRISMS_PF_BUILD_DEBUG} STREQUAL "ON")
        set_property(TARGET main_debug PROPERTY OUTPUT_NAME main-debug)
        deal_ii_setup_target(main_debug DEBUG)
        target_link_libraries(
            main_debug
            ${PRISMS_PF_CORE_DIR}/libprisms-pf-debug.a
        )

        if(${VTK_BUILT_SEPARATELY})
            include_directories(SYSTEM ${VTK_INCLUDE_DIR})
            target_link_libraries(main_debug ${VTK_NEW_LIBRARIES})
        endif()
        if(${PRISMS_PF_WITH_CALIPER})
            target_link_libraries(main_debug caliper)
        endif()
    endif()

    if(${PRISMS_PF_BUILD_RELEASE} STREQUAL "ON")
        set_property(TARGET main_release PROPERTY OUTPUT_NAME main)
        deal_ii_setup_target(main_release RELEASE)
        target_link_libraries(
            main_release
            ${PRISMS_PF_CORE_DIR}/libprisms-pf-release.a
        )

        if(${VTK_BUILT_SEPARATELY})
            include_directories(SYSTEM ${VTK_INCLUDE_DIR})
            target_link_libraries(main_release ${VTK_NEW_LIBRARIES})
        endif()
        if(${PRISMS_PF_WITH_CALIPER})
            target_link_libraries(main_release caliper)
        endif()
    endif()
endmacro()
