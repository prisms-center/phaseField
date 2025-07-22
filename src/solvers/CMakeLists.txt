# Manually specify files to be included
list(APPEND PRISMS_PF_SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/linear_solver_base.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/linear_solver_gmg.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/linear_solver_identity.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/nonexplicit_auxiliary_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/nonexplicit_base.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/nonexplicit_co_nonlinear_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/nonexplicit_linear_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/nonexplicit_self_nonlinear_solver.cc
)
set(PRISMS_PF_SOURCE_FILES ${PRISMS_PF_SOURCE_FILES} PARENT_SCOPE)