# Manually specify files to be included
target_sources(
  prisms_pf
  PRIVATE
    conditional_ostreams.cc
    constraint_manager.cc
    dof_manager.cc
    field_container.cc
    initial_conditions.cc
    nonuniform_dirichlet.cc
    problem.cc
    phase_field_tools.cc
    group_solution_handler.cc
    solution_indexer.cc
    solution_output.cc
    timer.cc
    triangulation_manager.cc
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS ${PROJECT_SOURCE_DIR}/include
      FILES
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/cell_marker_base.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/dependency_extents.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/field_container.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/invm_manager.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/phase_field_tools.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/solution_indexer.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/tee_stream.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/types.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/conditional_ostreams.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/dof_manager.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/grid_refiner_criterion.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/nonuniform_dirichlet.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/problem.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/solution_output.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/timer.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/constraint_manager.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/exceptions.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/group_solution_handler.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/parse_cmd_options.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/refinement_manager.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/solve_group.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/triangulation_manager.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/dependencies.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/field_attributes.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/initial_conditions.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/pde_operator_base.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/simulation_timer.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/system_wide.h
        ${PROJECT_SOURCE_DIR}/include/prismspf/core/type_enums.h
)

# Generated .inst files land in the build tree so .cc files
# that #include them need to find CMAKE_CURRENT_BINARY_DIR
target_include_directories(prisms_pf PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/..)

set(
  _inst_bases
  constraint_manager
  dof_manager
  field_container
  initial_conditions
  nonuniform_dirichlet
  problem
  phase_field_tools
  group_solution_handler
  solution_indexer
  solution_output
  triangulation_manager
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
