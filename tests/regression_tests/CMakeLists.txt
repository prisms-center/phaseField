##
# CMake for PRISMS-PF regression tests
##

message(STATUS "Setting up PRISMS-PF regression tests...")

# A list of the regression applications
set(REGRESSION_TEST_APPS
    allen_cahn_explicit
    allen_cahn_implicit
    cahn_hilliard_explicit
    cahn_hilliard_implicit
    cavity_flow
    fracture
    heat_equation_steady_state
    linear_solve_old_solution
    mechanics
    precipitate_explicit
    solution_blocks
)

# For each regression application install the source files, configure, and compile
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
        PATTERN "gold_output.txt"
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

# Run the regressions apps
find_package(Python COMPONENTS Interpreter Development REQUIRED)

if(NOT DEFINED N_THREADS)
 set(N_THREADS 1)
endif()

install(CODE "
execute_process(
  COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/run.py
    -j ${N_THREADS}
    -d ${CMAKE_INSTALL_PREFIX}/tests/regression_tests
    ${REGRESSION_TEST_APPS}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  RESULT_VARIABLE run_result
)
if(run_result)
  message(FATAL_ERROR \"Regression tests failed\")
endif()
")


message(STATUS "PRISMS-PF regression tests setup complete.")
