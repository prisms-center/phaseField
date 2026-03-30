# Manually specify files to be included
set(_sources)

set(
  _headers
  ${PROJECT_SOURCE_DIR}/include/prismspf/field_input/read_binary.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/field_input/read_field_base.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/field_input/read_field_factory.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/field_input/read_vtk.h
)

set(_inst_bases)

prisms_pf_target_sources(
  "${PRISMS_PF_TARGETS}"
  "${_headers}"
  "${_sources}"
  "${_inst_bases}"
)
