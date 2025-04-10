# Manually specify files to be included
list(APPEND PRISMS_PF_SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/conditional_ostreams.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/constraint_handler.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/dof_handler.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/initial_conditions.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/invm_handler.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/matrix_free_handler.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/matrix_free_operator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/nonuniform_dirichlet.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/pde_operator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/pde_problem.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solution_handler.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solution_output.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/timer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/triangulation_handler.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/variable_attribute_loader.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/variable_attributes.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/variable_container.cc
)
set(PRISMS_PF_SOURCE_FILES ${PRISMS_PF_SOURCE_FILES} PARENT_SCOPE)