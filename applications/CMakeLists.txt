# Manually list all the examples
set(
  _examples
  allen_cahn/explicit
  allen_cahn/implicit
  alloy_solidification
  blank
  cahn_hilliard_explicit
  coupled_allen_cahn_cahn_hilliard
  eshelby_inclusion
  nucleation
  laplace
)

# For each example install the source files, configure, and compile
foreach(_example ${_examples})
  # Install source files
  install(
    DIRECTORY
      ${_example}/
    DESTINATION applications/${_example}
    FILES_MATCHING
    PATTERN "CMakeFiles*"
      EXCLUDE
    PATTERN "solutions*"
      EXCLUDE
    PATTERN "*.cc"
    PATTERN "*.cpp"
    PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN "*.prm"
    PATTERN "CMakeLists.txt"
  )

  # Compile during installation
  # OUTPUT_QUIET on the configure step so we don't get extra verbosity.
  # CMAKE_INSTALL_MESSAGE NEVER so we don't print install and up-to-date
  # messages
  set(CMAKE_INSTALL_MESSAGE NEVER)
  install(
    CODE
      "
      execute_process(
        COMMAND ${CMAKE_COMMAND} -S . -B build
          -D PRISMS_PF_DIR=${CMAKE_BINARY_DIR}
        WORKING_DIRECTORY applications/${_example}
        RESULT_VARIABLE configure_result
        OUTPUT_QUIET
      )
      if(configure_result)
        message(FATAL_ERROR \"Failed to configure application: ${_example}\")
      else()
        execute_process(
          COMMAND ${CMAKE_COMMAND} --build build
          WORKING_DIRECTORY applications/${_example}
          RESULT_VARIABLE build_result
        )
        if(build_result)
          message(FATAL_ERROR \"Failed to build application: ${_example}\")
        endif()
      endif()
      "
  )
endforeach()
