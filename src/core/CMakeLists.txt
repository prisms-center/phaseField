# Manually specify files to be included
set(_src
    ${CMAKE_CURRENT_SOURCE_DIR}/conditional_ostreams.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/constraint_handler.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/dof_handler.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/initial_conditions.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/invm_handler.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/matrix_free_handler.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/matrix_free_operator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/nonuniform_dirichlet.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/pde_operator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/pde_problem.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/phase_field_tools.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/solution_handler.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/solution_output.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/solver_handler.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/timer.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/triangulation_handler.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/variable_attribute_loader.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/variable_attributes.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/variable_container.cc
)

set(_inst
    constraint_handler.inst.in
    dof_handler.inst.in
    initial_conditions.inst.in
    invm_handler.inst.in
    matrix_free_handler.inst.in
    matrix_free_operator.inst.in
    nonuniform_dirichlet.inst.in
    pde_operator.inst.in
    pde_problem.inst.in
    phase_field_tools.inst.in
    solution_handler.inst.in
    solution_output.inst.in
    solver_handler.inst.in
    triangulation_handler.inst.in
    variable_attribute_loader.inst.in
    variable_container.inst.in
)

file(
    GLOB _header
    CONFIGURE_DEPENDS
    ${CMAKE_SOURCE_DIR}/include/prismspf/core/*.h
)

define_library(prisms_pf_core OBJECT ${_src} ${_header} ${_inst})

# search the top level include directory
foreach(_build ${PRISMS_PF_BUILD_TYPES})
    string(TOLOWER ${_build} _build_lowercase)
    set(_target_name prisms_pf_core_${_build_lowercase})
    target_include_directories(
        ${_target_name}
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    )
    target_include_directories(
        ${_target_name}
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>
    )
endforeach()

expand_template_instantiations(prisms_pf_core "${_inst}")

collect_header_and_src_files("${_src}" "${_header}")
