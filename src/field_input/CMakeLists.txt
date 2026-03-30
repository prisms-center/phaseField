# Manually specify files to be included
target_sources(
  prisms_pf
  PRIVATE
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS ${PROJECT_SOURCE_DIR}/include
      FILES
        ${PROJECT_SOURCE_DIR}/include/prismspf/field_input/read_binary.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/field_input/read_field_base.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/field_input/read_field_factory.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/field_input/read_vtk.h
)
