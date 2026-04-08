# Manually specify files to be included
set(_sources integrator.cc)

set(
  _headers
  ${PROJECT_SOURCE_DIR}/include/prismspf/utilities/integrator.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/utilities/periodic_distance.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/utilities/utilities.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/utilities/vectorized_operations.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/utilities/crystal_symmetry.h
)

set(_inst_bases integrator)

prisms_pf_target_sources(
  "${PRISMS_PF_TARGETS}"
  "${_headers}"
  "${_sources}"
  "${_inst_bases}"
)
