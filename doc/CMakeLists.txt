#
# CMake for creating the documentation. There's two parts.
#   1. Scrape the parameters to build documentation
#   2. Build doxygen documentation
#

# First, we must scrape the parameters to build their documentation
add_subdirectory(scripts)
add_custom_target(generate_parameter_docs
    COMMAND scrape_parameters_exe
    COMMENT "Generating parameter documentation"
)
add_dependencies(generate_parameter_docs scrape_parameters_exe)


find_package(Doxygen REQUIRED)

# Path of Doxyfile and its output
set(PRISMS_PF_DOXYGEN_INPUT ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile)
set(PRISMS_PF_DOXYGEN_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen/Doxyfile.out)

# The directory where the documentation will be generated
set(DOXYGEN_HTML_OUTPUT "${PRISMS_PF_VERSION}")
set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Copy input to output
configure_file(${PRISMS_PF_DOXYGEN_INPUT} ${PRISMS_PF_DOXYGEN_OUTPUT} @ONLY)

# Add target to generate documentation
doxygen_add_docs(
  doc
  COMMENT "Generating Doxygen documentation"
  CONFIG_FILE ${PRISMS_PF_DOXYGEN_OUTPUT}
  ALL
)
add_dependencies(doc generate_parameter_docs)

message(
  STATUS
  "Doxygen documentation will be generated using:\n"
  "     ${PRISMS_PF_DOXYGEN_INPUT}\n"
  "   and output to:\n"
  "     ${DOXYGEN_OUTPUT_DIRECTORY}/${DOXYGEN_HTML_OUTPUT}\n"
)
