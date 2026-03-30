# Manually specify files to be included
target_sources(
  prisms_pf
  PRIVATE
    input_file_reader.cc
    user_input_parameters.cc
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS ${PROJECT_SOURCE_DIR}/include
      FILES
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

# Generated .inst files land in the build tree so .cc files
# that #include them need to find CMAKE_CURRENT_BINARY_DIR
target_include_directories(prisms_pf PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/..)

set(_inst_bases user_input_parameters)

# Derive the .inst.in list from the base names
foreach(_base ${_inst_bases})
  list(APPEND _inst "${_base}.inst.in")
endforeach()

expand_template_instantiations(prisms_pf "${_inst}")

# Let CMake know that .cc files can't compile until the .inst
# files exist
foreach(_base ${_inst_bases})
  set_source_files_properties(
    ${_base}.cc
    PROPERTIES
      OBJECT_DEPENDS
        "${CMAKE_CURRENT_BINARY_DIR}/${_base}.inst"
  )
endforeach()
