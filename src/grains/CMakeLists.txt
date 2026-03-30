# Manually specify files to be included
set(_sources)

set(_headers ${PROJECT_SOURCE_DIR}/include/prismspf/grains/grains.h)

set(_inst_bases)

prisms_pf_target_sources(
  "${PRISMS_PF_TARGETS}"
  "${_headers}"
  "${_sources}"
  "${_inst_bases}"
)
