# Manually specify files to be included
set(_src
    ${CMAKE_CURRENT_SOURCE_DIR}/integrator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/element_volume.cc
)

set(_inst integrator.inst.in element_volume.inst.in)

foreach(_build ${PRISMS_PF_BUILD_TYPES})
    string(TOLOWER ${_build} _build_lowercase)
    add_library(prisms_pf_utilities_${_build_lowercase} STATIC ${_src})
    target_include_directories(
        prisms_pf_utilities_${_build_lowercase}
        PRIVATE ${CMAKE_BINARY_DIR}/include ${CMAKE_BINARY_DIR}/src
    )
    if(${VTK_BUILT_SEPARATELY})
        include_directories(SYSTEM ${VTK_INCLUDE_DIR})
        target_link_libraries(
            prisms_pf_utilities_${_build_lowercase}
            ${VTK_NEW_LIBRARIES}
        )
    endif()
    if(${PRISMS_PF_WITH_CALIPER})
        target_link_libraries(prisms_pf_utilities_${_build_lowercase} caliper)
    endif()
    deal_ii_setup_target(prisms_pf_utilities_${_build_lowercase} DEBUG)
endforeach()

expand_template_instantiations(prisms_pf_utilities "${_inst}")

list(APPEND PRISMS_PF_SOURCE_FILES ${_src})

set(PRISMS_PF_SOURCE_FILES ${PRISMS_PF_SOURCE_FILES} PARENT_SCOPE)
