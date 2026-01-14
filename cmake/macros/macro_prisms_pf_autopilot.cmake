#
# Autopilot function that will set up an application
#

macro(prisms_pf_autopilot PRISMS_PF_CORE_DIR)
    # Add the script files
    add_subdirectory(
        "${PRISMS_PF_CORE_DIR}/lib/cmake/prisms_pf/scripts"
        "${CMAKE_BINARY_DIR}/scripts_build"
    )

    # Enable compile commands export
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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

    # Create the exectuables
    foreach(_build ${APPLICATION_BUILD_TYPES})
        string(TOLOWER ${_build} _build_lowercase)
        string(TOUPPER ${_build} _build_uppercase)

        set(_target main-${_build_lowercase})

        add_executable(${_target} ${TARGET_SRC})

        if(TARGET prisms_pf::prisms_pf_${_build_lowercase})
            target_link_libraries(
                ${_target}
                PRIVATE prisms_pf::prisms_pf_${_build_lowercase}
            )
            message(
                STATUS
                "Successfully linked ${_target} to prisms_pf_${_build_lowercase}"
            )
        else()
            message(
                FATAL_ERROR
                "Target prisms_pf_${_build_lowercase} does not exist! "
                "Available build types: ${PRISMS_PF_BUILD_TYPES}"
            )
        endif()

        target_include_directories(
            ${_target}
            PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
        )

        target_compile_options(
            ${_target}
            PRIVATE
                $<$<COMPILE_LANGUAGE:CXX>:
                ${PRISMS_PF_WARNING_FLAGS}
                ${PRISMS_PF_CXX_FLAGS}
                ${PRISMS_PF_CXX_FLAGS_${_build_uppercase}}>
        )
        target_link_options(
            ${_target}
            PRIVATE
                $<$<LINK_LANGUAGE:CXX>:
                ${PRISMS_PF_LINKER_FLAGS}
                ${PRISMS_PF_LINKER_FLAGS_${_build_uppercase}}>
        )

        set(BUILD_OVERRIDE "")
        expand_template_instantiations(${_target} "custom_pde.inst.in")
    endforeach()
endmacro()
