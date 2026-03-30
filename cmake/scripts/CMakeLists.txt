# Add the script to an executable
add_executable(expand_template_instantiations_exe expand_template_instantiations.cc)

target_compile_features(expand_template_instantiations_exe PRIVATE cxx_std_20)

set_target_properties(
  expand_template_instantiations_exe
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
      "${CMAKE_BINARY_DIR}/cmake/scripts"
    CXX_EXTENSIONS
      OFF
)
