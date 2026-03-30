# Manually specify files to be included
set(
  _sources
  input_file_reader.cc
  user_input_parameters.cc
)

set(
  _headers
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/checkpoint_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/load_initial_condition_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/spatial_discretization.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/constraint_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/miscellaneous_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/temporal_discretization.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/input_file_reader.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/nonlinear_solve_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/user_constants.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/io_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/nucleation_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/user_input_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/linear_solve_parameters.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/user_inputs/output_parameters.h
)

set(_inst_bases user_input_parameters)

prisms_pf_target_sources(
  "${PRISMS_PF_TARGETS}"
  "${_headers}"
  "${_sources}"
  "${_inst_bases}"
)
