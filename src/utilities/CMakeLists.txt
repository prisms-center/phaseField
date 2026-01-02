# Manually specify files to be included
set(
  _src
  ${CMAKE_CURRENT_SOURCE_DIR}/integrator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/element_volume.cc
)

set(
  _inst
  integrator.inst.in
  element_volume.inst.in
)

file(GLOB _header CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/include/prismspf/utilities/*.h)

define_library(
  prisms_pf_utilities
  OBJECT
  ${_src}
  ${_header}
  ${_inst}
)

expand_template_instantiations(prisms_pf_utilities "${_inst}")

collect_header_and_src_files("${_src}" "${_header}")
