# Manually list all the examples
set(
  _examples
  allen_cahn/explicit
  allen_cahn/implicit
  alloy_solidification
  blank
  cahn_hilliard/explicit
  cahn_hilliard/implicit
  coupled_allen_cahn_cahn_hilliard
  mechanics/boundary_value_problem
  mechanics/eshelby_inclusion
  nucleation
  laplace
  pfhub_benchmarks/1/a
  pfhub_benchmarks/2/b
  pfhub_benchmarks/6/a
)

# Create a target for each example and feed them into a big target
set(_example_targets)

foreach(_example IN LISTS _examples)
  string(REPLACE "/" "_" _example_name ${_example})

  add_custom_target(
    ${_example_name}
    COMMAND
      ${CMAKE_COMMAND} --log-level=NOTICE -S ${CMAKE_SOURCE_DIR}/applications/${_example}
      -B ${CMAKE_BINARY_DIR}/applications/${_example}
      -DCMAKE_PREFIX_PATH=${CMAKE_INSTALL_PREFIX} -DPRISMS_PF_DIR=${CMAKE_BINARY_DIR}
    COMMAND
      ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}/applications/${_example}
    DEPENDS
      ${PRISMS_PF_TARGETS}
  )

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

  install(
    PROGRAMS
      ${CMAKE_BINARY_DIR}/applications/${_example}/main
    DESTINATION applications/${_example}
  )

  list(APPEND _example_targets ${_example_name})
endforeach()

add_custom_target(build_examples ALL DEPENDS ${_example_targets})
