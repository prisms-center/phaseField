# Add the script to an executable
add_executable(scrape_parameters_exe scrape_parameters.cc)

target_compile_features(scrape_parameters_exe PRIVATE cxx_std_20)

target_link_libraries(scrape_parameters_exe PRIVATE prisms_pf)

set_target_properties(
  scrape_parameters_exe
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}"
    CXX_EXTENSIONS
      OFF
)
