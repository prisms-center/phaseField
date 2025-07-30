# Set the c++20 standard
set(CMAKE_CXX_STANDARD 20)

# Add the script to an executable
add_executable(expand_template_instantiations_exe
  ${CMAKE_CURRENT_SOURCE_DIR}/expand_template_instantiations.cc
)

set_target_properties(expand_template_instantiations_exe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cmake/scripts"
)
