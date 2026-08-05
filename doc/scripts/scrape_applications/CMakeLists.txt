# Add the script to an executable
add_executable(scrape_applications_exe scrape_applications.cc)

target_compile_features(scrape_applications_exe PRIVATE cxx_std_20)

target_link_libraries(scrape_applications_exe PRIVATE ${PRISMS_PF_TARGETS})

set_target_properties(
  scrape_applications_exe
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}"
    CXX_EXTENSIONS
      OFF
)
