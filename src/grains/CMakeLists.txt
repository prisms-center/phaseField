# Manually specify files to be included
target_sources(
  prisms_pf
  PRIVATE
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS ${PROJECT_SOURCE_DIR}/include
      FILES ${PROJECT_SOURCE_DIR}/include/prismspf/grains/grains.h
)
