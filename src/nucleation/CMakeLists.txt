# Manually specify files to be included
target_sources(
  prisms_pf
  PRIVATE
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS ${PROJECT_SOURCE_DIR}/include
      FILES
        ${PROJECT_SOURCE_DIR}/include/prismspf/nucleation/nucleation_manager.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/nucleation/nucleus.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/nucleation/nucleus_refinement_function.h
)
