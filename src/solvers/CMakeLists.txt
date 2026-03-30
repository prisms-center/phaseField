# Manually specify files to be included
set(
  _sources
  mf_operator.cc
  solver_base.cc
  solvers.cc
)

set(
  _headers
  ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/constant_solver.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/explicit_solver.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/linear_helper.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/linear_solver.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/mf_operator.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/mg_solver.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/newton_solver.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/solve_context.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/solver_base.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/solvers.h
)

set(
  _inst_bases
  mf_operator
  solver_base
  solvers
)

prisms_pf_target_sources(
  "${PRISMS_PF_TARGETS}"
  "${_headers}"
  "${_sources}"
  "${_inst_bases}"
)
