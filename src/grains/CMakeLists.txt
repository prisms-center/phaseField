# Manually specify files to be included
list(APPEND PRISMS_PF_SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/FloodFiller.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/OrderParameterRemapper.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/reassignGrains.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/SimplifiedGrainRepresentation.cc
)
set(PRISMS_PF_SOURCE_FILES ${PRISMS_PF_SOURCE_FILES} PARENT_SCOPE)