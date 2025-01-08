# Manually specify files to be included
list(APPEND PRISMS_PF_SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/user_inputs/input_file_reader.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/user_inputs/user_input_parameters.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/conditional_ostreams.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/variable_attribute_loader.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/variable_attributes.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/variable_container.cc
)
set(PRISMS_PF_SOURCE_FILES ${PRISMS_PF_SOURCE_FILES} PARENT_SCOPE)