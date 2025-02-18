# Manually specify files to be included
list(APPEND PRISMS_PF_SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/input_file_reader.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/user_input_parameters.cc
)
set(PRISMS_PF_SOURCE_FILES ${PRISMS_PF_SOURCE_FILES} PARENT_SCOPE)