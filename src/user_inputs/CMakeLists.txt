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

# search the top level include directory
foreach(_build ${PRISMS_PF_BUILD_TYPES})
    string(TOLOWER ${_build} _build_lowercase)
    set(_target_name prisms_pf_user_inputs_${_build_lowercase})
    target_include_directories(
        ${_target_name}
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    )
    target_include_directories(
        ${_target_name}
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>
    )
endforeach()

expand_template_instantiations(prisms_pf_user_inputs "${_inst}")

collect_header_and_src_files("${_src}" "${_header}")
