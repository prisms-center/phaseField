# Manually specify files to be included
set(
  _sources
  flood_filler.cc
  grains.cc
  grain_reassignment.cc
)

set(
  _headers
  ${PROJECT_SOURCE_DIR}/include/prismspf/grains/flood_filler.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/grains/grains.h
  ${PROJECT_SOURCE_DIR}/include/prismspf/grains/grain_reassignment.h
)

set(
  _inst_bases
  flood_filler
  grains
  grain_reassignment
)

prisms_pf_target_sources(
  "${PRISMS_PF_TARGETS}"
  "${_headers}"
  "${_sources}"
  "${_inst_bases}"
)
