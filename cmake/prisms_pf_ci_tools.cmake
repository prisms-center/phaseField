#
# CI tools
#

if(PRISMS_PF_CLANG_TIDY)
  # NOTE: I've tested clang-tidy-18 and found that it fails. clang-tidy-20 works.
  find_program(
    CLANG_TIDY_TOOL
    NAMES
      clang-tidy
      clang-tidy-20
      clang-tidy-19
  )
  execute_process(
    COMMAND
      ${CLANG_TIDY_TOOL} --version
    OUTPUT_VARIABLE CLANG_TIDY_TOOL_VERSION
    ERROR_VARIABLE CLANG_TIDY_TOOL_VERSION
  )
  string(
    REGEX MATCH
    "[0-9]+(\\.[0-9]+)+"
    CLANG_TIDY_TOOL_VERSION
    "${CLANG_TIDY_TOOL_VERSION}"
  )
  message(STATUS "clang-tidy ${CLANG_TIDY_TOOL_VERSION} (${CLANG_TIDY_TOOL})")
  # NOTE: To avoid linting on the vendored source files we set CXX_CLANG_TIDY
  # as a property of our targets
endif()
