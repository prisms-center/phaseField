# Manually specify files to be included
set(_src
    ${CMAKE_CURRENT_SOURCE_DIR}/input_file_reader.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/user_input_parameters.cc
)

set(_inst user_input_parameters.inst.in)

file(
    GLOB _header
    CONFIGURE_DEPENDS
    ${CMAKE_SOURCE_DIR}/include/prismspf/user_inputs/*.h
)

define_library(prisms_pf_user_inputs OBJECT ${_src} ${_header} ${_inst})

expand_template_instantiations(prisms_pf_user_inputs "${_inst}")

collect_header_and_src_files("${_src}" "${_header}")
