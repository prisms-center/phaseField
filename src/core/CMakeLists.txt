# Manually specify files to be included
list(APPEND PRISMS_PF_SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/boundaryConditions.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/markBoundaries.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions/vectorBCFunction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/initial_conditions/initialConditions.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/postprocessing/computeIntegral.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/postprocessing/postprocessor.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/refinement/AdaptiveRefinement.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/computeLHS.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/computeRHS.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/setNonlinearEqInitialGuess.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/solve.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/solveIncrement.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/SolverParameters.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/buildFields.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/checkpoint.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/init.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/inputFileReader.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/invM.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/matrixFreePDE.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/outputResults.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/reinit.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/userInputParameters.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/variableAttributeLoader.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/variableAttributes.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/variableContainer.cc
)
set(PRISMS_PF_SOURCE_FILES ${PRISMS_PF_SOURCE_FILES} PARENT_SCOPE)