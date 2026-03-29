# Manually specify files to be included
target_sources(
  prisms_pf
  PRIVATE
    integrator.cc
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS ${PROJECT_SOURCE_DIR}/include
      FILES
        ${PROJECT_SOURCE_DIR}/include/prismspf/utilities/integrator.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/utilities/periodic_distance.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/utilities/utilities.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/utilities/vectorized_operations.h
)

# Generated .inst files land in the build tree so .cc files
# that #include them need to find CMAKE_CURRENT_BINARY_DIR
target_include_directories(prisms_pf PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/..)

set(_inst_bases integrator)

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
