# Manually specify files to be included
target_sources(
  prisms_pf
  PRIVATE
    mf_operator.cc
    solver_base.cc
    solvers.cc
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS ${PROJECT_SOURCE_DIR}/include
      FILES
        ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/constant_solver.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/explicit_solver.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/linear_helper.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/linear_solver.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/mf_operator.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/mg_solver.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/newton_solver.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/solve_context.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/solver_base.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/solvers/solvers.h
)

# Generated .inst files land in the build tree so .cc files
# that #include them need to find CMAKE_CURRENT_BINARY_DIR
target_include_directories(prisms_pf PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/..)

set(
  _inst_bases
  mf_operator
  solver_base
  solvers
)

# Derive the .inst.in list from the base names
foreach(_base ${_inst_bases})
  list(APPEND _inst "${_base}.inst.in")
endforeach()

expand_template_instantiations(prisms_pf "${_inst}")

# Let CMake know that .cc files can't compile until the .inst
# files exist
foreach(_base ${_inst_bases})
  set_source_files_properties(
    ${_base}.cc
    PROPERTIES
      OBJECT_DEPENDS
        "${CMAKE_CURRENT_BINARY_DIR}/${_base}.inst"
  )
endforeach()
