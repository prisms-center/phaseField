# Manually specify files to be included
set(_src
    ${CMAKE_CURRENT_SOURCE_DIR}/group_solver_base.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/group_solvers.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/mf_operator.cc
)

set(_inst
    group_solver_base.inst.in
    group_solvers.inst.in
    mf_operator.inst.in
)

file(GLOB _header CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/include/prismspf/solvers/*.h)

define_library(
  prisms_pf_solvers
  OBJECT
  ${_src}
  ${_header}
  ${_inst}
)

expand_template_instantiations(prisms_pf_solvers "${_inst}")

collect_header_and_src_files("${_src}" "${_header}")
