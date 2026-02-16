# Manually specify files to be included
set(_src
    ${CMAKE_CURRENT_SOURCE_DIR}/conditional_ostreams.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/constraint_manager.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/dof_manager.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/initial_conditions.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/nonuniform_dirichlet.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/pde_operator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/problem.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/phase_field_tools.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/solution_indexer.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/solution_output.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/timer.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/triangulation_manager.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/variable_container.cc
)

set(_inst
    constraint_manager.inst.in
    dof_manager.inst.in
    initial_conditions.inst.in
    nonuniform_dirichlet.inst.in
    pde_operator.inst.in
    problem.inst.in
    phase_field_tools.inst.in
    solution_indexer.inst.in
    solution_output.inst.in
    triangulation_manager.inst.in
)

file(GLOB _header CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/include/prismspf/core/*.h)

define_library(
  prisms_pf_core
  OBJECT
  ${_src}
  ${_header}
  ${_inst}
)

expand_template_instantiations(prisms_pf_core "${_inst}")

collect_header_and_src_files("${_src}" "${_header}")
