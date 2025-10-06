##
# CMake for PRISMS-PF regression tests
##

message(STATUS "Setting up PRISMS-PF regression tests...")

set(REGRESSION_TEST_APPS
    allen_cahn_explicit
    allen_cahn_implicit
    allen_cahn_rk4
    blank
    cahn_hilliard_explicit
    cahn_hilliard_implicit
    cavity_flow
    file_read_binary
    file_read_vtk
    fracture
    heat_equation_steady_state
    linear_solve_old_solution
    mechanics
    precipitate_explicit
    solution_blocks
    taylor_green_vortex
)

foreach(test_app ${REGRESSION_TEST_APPS})
    # Install the source files
    install(
        DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${test_app}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/tests/regression_tests/${test_app}
        FILES_MATCHING
        PATTERN "CMakeFiles*" EXCLUDE
        PATTERN "*.in"
        PATTERN "*.cc"
        PATTERN "*.h"
        PATTERN "*.prm"
        PATTERN "CMakeLists.txt"
    )

    # Compile during installation
    install(
        CODE
            "
    message(STATUS \"Building regression test: ${test_app}\")
    execute_process(
      COMMAND ${CMAKE_COMMAND} 
        -DCMAKE_BUILD_TYPE=Release
        -S ${CMAKE_INSTALL_PREFIX}/tests/regression_tests/${test_app}/
      WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/tests/regression_tests/${test_app}
      RESULT_VARIABLE configure_result
    )
    if(configure_result)
      message(WARNING \"Failed to configure regression test: ${test_app}\")
    else()
        execute_process(
        COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/tests/regression_tests/${test_app}
        RESULT_VARIABLE compile_result
        )
        if(compile_result)
            message(WARNING \"Failed to compile regression test: ${test_app}\")
        endif()
    endif()
  "
    )
endforeach()

message(STATUS "PRISMS-PF regression tests setup complete.")
