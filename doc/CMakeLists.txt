#
# CMake for creating the documentation. There's two parts.
#   1. Scrape the parameters to build documentation
#   2. Build doxygen documentation
#

add_subdirectory(scripts)

# First, we must scrape the parameters to build their documentation
add_custom_target(
  generate_parameter_docs
  COMMAND
    scrape_parameters_exe ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/parameters/parameters.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/parameters/parameters.h
  COMMENT "Generating parameter documentation"
)
add_dependencies(generate_parameter_docs scrape_parameters_exe)

# Next, scrape the applications to build their documentation
add_custom_target(
  generate_application_docs
  COMMAND
    scrape_applications_exe
    ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/applications/applications.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/applications/applications.h
    ${CMAKE_SOURCE_DIR}/applications
  COMMENT "Generating parameter documentation"
)
add_dependencies(generate_application_docs scrape_applications_exe)

# Next, doxygen documentation with the parameter docs dependency
find_package(Doxygen REQUIRED)

# Path of Doxyfile and its output
set(PRISMS_PF_DOXYGEN_INPUT ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile)
set(PRISMS_PF_DOXYGEN_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen/Doxyfile.out)

# The directory where the documentation will be generated
set(DOXYGEN_HTML_OUTPUT "${PRISMS_PF_GIT_BRANCH}")
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
add_dependencies(
  doc
  generate_parameter_docs
  generate_application_docs
)

message(
  STATUS
  "Doxygen documentation will be generated using:\n"
  "     ${PRISMS_PF_DOXYGEN_INPUT}\n"
  "   and output to:\n"
  "     ${DOXYGEN_OUTPUT_DIRECTORY}/${DOXYGEN_HTML_OUTPUT}\n"
)
