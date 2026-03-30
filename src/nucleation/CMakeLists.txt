# Manually specify files to be included
set(_sources)

set(
  _headers
  ${PROJECT_SOURCE_DIR}/include/prismspf/nucleation/nucleation_manager.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/nucleation/nucleus.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/nucleation/nucleus_refinement_function.h
)

set(_inst_bases)

prisms_pf_target_sources(
  "${PRISMS_PF_TARGETS}"
  "${_headers}"
  "${_sources}"
  "${_inst_bases}"
)
