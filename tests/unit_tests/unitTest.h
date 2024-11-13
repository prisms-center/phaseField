#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/timer.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

// define data type
template <int dim>
void
computeStress(const dealii::Table<2, double>       &CIJ,
              const dealii::VectorizedArray<double> ux[][dim],
              const dealii::VectorizedArray<double> R[][dim]);

#include "../../include/matrixFreePDE.h"
#include "../../include/parallelNucleationList.h"
#include "../../src/FloodFiller/FloodFiller.cc"
#include "../../src/OrderParameterRemapper/OrderParameterRemapper.cc"
#include "../../src/SimplifiedGrainRepresentation/SimplifiedGrainRepresentation.cc"
#include "../../src/inputFileReader/inputFileReader.cc"
#include "../../src/matrixfree/AdaptiveRefinement.cc"
#include "../../src/matrixfree/boundaryConditions.cc"
#include "../../src/matrixfree/checkpoint.cc"
#include "../../src/matrixfree/computeIntegral.cc"
#include "../../src/matrixfree/computeLHS.cc"
#include "../../src/matrixfree/computeRHS.cc"
#include "../../src/matrixfree/init.cc"
#include "../../src/matrixfree/initForTests.cc"
#include "../../src/matrixfree/initialConditions.cc"
#include "../../src/matrixfree/invM.cc"
#include "../../src/matrixfree/markBoundaries.cc"
#include "../../src/matrixfree/matrixFreePDE.cc"
#include "../../src/matrixfree/nucleation.cc"
#include "../../src/matrixfree/outputResults.cc"
#include "../../src/matrixfree/postprocessor.cc"
#include "../../src/matrixfree/reassignGrains.cc"
#include "../../src/matrixfree/reinit.cc"
#include "../../src/matrixfree/setNonlinearEqInitialGuess.cc"
#include "../../src/matrixfree/solve.cc"
#include "../../src/matrixfree/solveIncrement.cc"
#include "../../src/matrixfree/utilities.cc"
#include "../../src/models/mechanics/computeStress.h"
#include "../../src/parallelNucleationList/parallelNucleationList.cc"
#include "../../src/userInputParameters/loadVariableAttributes.cc"
#include "../../src/userInputParameters/load_BC_list.cc"
#include "../../src/userInputParameters/load_user_constants.cc"
#include "../../src/userInputParameters/setTimeStepList.cc"
#include "../../src/variableAttributeLoader/variableAttributeLoader.cc"
#include "../../src/variableContainer/variableContainer.cc"

template <int dim, typename T>
class unitTest
{
public:
  bool
  test_computeInvM(int argc, char **argv, userInputParameters<dim>);
  bool
  test_outputResults(int argc, char **argv, userInputParameters<dim> userInputs);
  bool
  test_computeStress();
  void
  assignCIJSize(
    dealii::VectorizedArray<double> CIJ[2 * dim - 1 + dim / 3][2 * dim - 1 + dim / 3]);
  void
  assignCIJSize(dealii::Table<2, double> &CIJ);
  bool
  test_setRigidBodyModeConstraints(std::vector<int>, userInputParameters<dim> userInputs);
  bool
  test_parse_line();
  bool
  test_get_subsection_entry_list();
  bool
  test_get_entry_name_ending_list();
  bool
  test_load_BC_list();
  bool
  test_setOutputTimeSteps();
  bool
  test_NonlinearSolverParameters();
  bool
  test_LinearSolverParameters();
  bool
  test_EquationDependencyParser_variables_and_residuals_needed();
  bool
  test_EquationDependencyParser_nonlinear();
  bool
  test_EquationDependencyParser_postprocessing();
  bool
  test_FloodFiller();
  bool
  test_SimplifiedGrainRepresentation();
  bool
  test_SimplifiedGrainManipulator_transferGrainIds();
  bool
  test_SimplifiedGrainManipulator_reassignGrains();
  bool
  test_OrderParameterRemapper();
};

#include "test_FloodFiller.h"
#include "test_LinearSolverParameters.h"
#include "test_NonlinearSolverParameters.h"
#include "test_OrderParameterRemapper.h"
#include "test_SimplifiedGrainManipulator.h"
#include "test_SimplifiedGrainRepresentation.h"
#include "test_computeStress.h"
#include "test_get_entry_name_ending_list.h"
#include "test_get_subsection_entry_list.h"
#include "test_invM.h"
#include "test_load_BC_list.h"
#include "test_outputResults.h"
#include "test_parse_line.h"
#include "test_setOutputTimeSteps.h"
#include "variableAttributeLoader_test.cc"

#include "../../include/SolverParameters.h"
#include "../../src/SolverParameters/SolverParameters.cc"
