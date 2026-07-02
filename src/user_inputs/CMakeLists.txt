# Manually specify files to be included
set(
  _sources
  constraint_parameters.cc
  io_parameters.cc
  miscellaneous_parameters.cc
  nucleation_parameters.cc
  solve_parameters.cc
  spatial_discretization.cc
  temporal_discretization.cc
  user_input_parameters.cc
)

set(
  _headers
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/constraint_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/io_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/miscellaneous_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/nucleation_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/parameter_base.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/solve_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/spatial_discretization.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/temporal_discretization.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/user_constants.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/user_input_parameters.h
)

set(
  _inst_bases
  constraint_parameters
  nucleation_parameters
  spatial_discretization
  user_input_parameters
)

prisms_pf_target_sources(
  "${PRISMS_PF_TARGETS}"
  "${_headers}"
  "${_sources}"
  "${_inst_bases}"
)
