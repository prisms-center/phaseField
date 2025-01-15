# Manually specify files to be included
list(APPEND PRISMS_PF_SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/nucleation.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/parallelNucleationList.cc
)
set(PRISMS_PF_SOURCE_FILES ${PRISMS_PF_SOURCE_FILES} PARENT_SCOPE)