# Manually specify files to be included
set(
  _src
  ${CMAKE_CURRENT_SOURCE_DIR}/concurrent_constant_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/concurrent_explicit_postprocess_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/concurrent_explicit_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/concurrent_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/linear_solver_base.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/linear_solver_gmg.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/linear_solver_identity.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/sequential_auxiliary_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/sequential_co_nonlinear_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/sequential_linear_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/sequential_self_nonlinear_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/sequential_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solver_base.cc
)

set(
  _inst
  concurrent_constant_solver.inst.in
  concurrent_explicit_postprocess_solver.inst.in
  concurrent_explicit_solver.inst.in
  concurrent_solver.inst.in
  linear_solver_base.inst.in
  linear_solver_gmg.inst.in
  linear_solver_identity.inst.in
  sequential_auxiliary_solver.inst.in
  sequential_co_nonlinear_solver.inst.in
  sequential_linear_solver.inst.in
  sequential_self_nonlinear_solver.inst.in
  sequential_solver.inst.in
  solver_base.inst.in
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
