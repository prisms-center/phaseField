# Add the script to an executable
add_executable(scrape_parameters_exe scrape_parameters.cc)

prisms_pf_configure_target(scrape_parameters_exe "Debug")
target_include_directories(scrape_paramteres_exe PUBLIC
    ${PROJECT_SOURCE_DIR}/include/prismspf
)

set_target_properties(
  scrape_parameters_exe
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
		"${CMAKE_CURRENT_BINARY_DIR}"
    CXX_EXTENSIONS
      OFF
)
