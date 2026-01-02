##
# CMake for PRISMS-PF applications
##

message(STATUS "Setting up PRISMS-PF applications...")

# A list of the applications
set(
  APPLICATIONS
  allen_cahn_explicit
  alloy_solidification
  blank
  cahn_hilliard_explicit
  coupled_allen_cahn_cahn_hilliard
  eshelby_inclusion
  MgNd_precipitate_single_Bppp
  nucleation
  precipitate_evolution
  spinodal_decomposition
)

# For each regression application install the source files, configure, and compile
foreach(app ${APPLICATIONS})
  # Install the source files
  install(
    DIRECTORY
      ${CMAKE_CURRENT_SOURCE_DIR}/${app}/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/applications/${app}
    FILES_MATCHING
    PATTERN "CMakeFiles*"
      EXCLUDE
    PATTERN "*.in"
    PATTERN "*.cc"
    PATTERN "*.h"
    PATTERN "*.prm"
    PATTERN "CMakeLists.txt"
    PATTERN "*.md"
  )

  # Compile during installation
  install(
    CODE
      "
    message(STATUS \"Building regression test: ${app}\")
    execute_process(
      COMMAND ${CMAKE_COMMAND}
        -S ${CMAKE_INSTALL_PREFIX}/applications/${app}/
      WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/applications/${app}
      RESULT_VARIABLE configure_result
    )
    if(configure_result)
      message(WARNING \"Failed to configure application: ${app}\")
    else()
        execute_process(
        COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/applications/${app}
        RESULT_VARIABLE compile_result
        )
        if(compile_result)
            message(WARNING \"Failed to compile application: ${app}\")
        endif()
    endif()
  "
  )
endforeach()

message(STATUS "PRISMS-PF applications setup complete.")
