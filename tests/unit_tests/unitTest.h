#include "../../include/dealIIheaders.h"
#include <iostream>

//define data type
template <int dim>
void computeStress(const dealii::Table<2, double>& CIJ, const dealii::VectorizedArray<double> ux[][dim], const dealii::VectorizedArray<double> R[][dim]);

#include "../../include/matrixFreePDE.h"
#include "../../include/parallelNucleationList.h"

#include "../../src/matrixfree/matrixFreePDE.cc"
#include "../../src/matrixfree/init.cc"
#include "../../src/matrixfree/reinit.cc"
#include "../../src/matrixfree/initForTests.cc"
#include "../../src/matrixfree/refine.cc"
#include "../../src/matrixfree/invM.cc"
#include "../../src/matrixfree/computeLHS.cc"
#include "../../src/matrixfree/computeRHS.cc"
#include "../../src/matrixfree/solve.cc"
#include "../../src/matrixfree/solveIncrement.cc"
#include "../../src/matrixfree/outputResults.cc"
#include "../../src/matrixfree/markBoundaries.cc"
#include "../../src/matrixfree/boundaryConditions.cc"
#include "../../src/matrixfree/initialConditions.cc"
#include "../../src/matrixfree/utilities.cc"
#include "../../src/matrixfree/computeIntegral.cc"
#include "../../src/matrixfree/nucleation.cc"
#include "../../src/matrixfree/checkpoint.cc"

#include "../../src/matrixfree/setNonlinearEqInitialGuess.cc"

#include "../../src/inputFileReader/inputFileReader.cc"
#include "../../src/parallelNucleationList/parallelNucleationList.cc"

#include "../../src/models/mechanics/computeStress.h"
#include "../../src/matrixfree/postprocessor.cc"

#include "../../src/utilities/sortIndexEntryPairList.cc"

#include "../../src/variableContainer/variableContainer.cc"

#include "../../src/userInputParameters/load_BC_list.cc"
#include "../../src/userInputParameters/load_user_constants.cc"
#include "../../src/userInputParameters/loadVariableAttributes.cc"
#include "../../src/userInputParameters/setTimeStepList.cc"

#include "../../src/variableAttributeLoader/variableAttributeLoader.cc"

template <int dim, typename T>
class unitTest
{
	public:
	bool test_computeInvM(int argc, char **argv, userInputParameters<dim>);
	bool test_outputResults(int argc, char **argv, userInputParameters<dim> userInputs);
	bool test_computeStress();
	void assignCIJSize(dealii::VectorizedArray<double> CIJ[2*dim-1+dim/3][2*dim-1+dim/3]);
	void assignCIJSize(dealii::Table<2, double> &CIJ);
	bool test_setRigidBodyModeConstraints(std::vector<int>, userInputParameters<dim> userInputs);
	bool test_parse_line();
	bool test_get_subsection_entry_list();
	bool test_get_entry_name_ending_list();
	bool test_load_BC_list();
	bool test_setOutputTimeSteps();
    bool test_NonlinearSolverParameters();
    bool test_LinearSolverParameters();
    bool test_EquationDependencyParser_variables_and_residuals_needed();
    bool test_EquationDependencyParser_nonlinear();
    bool test_EquationDependencyParser_postprocessing();
};

#include "variableAttributeLoader_test.cc"
#include "test_invM.h"
#include "test_outputResults.h"
#include "test_computeStress.h"
#include "test_setOutputTimeSteps.h"
#include "test_setRigidBodyModeConstraints.h"
#include "test_parse_line.h"
#include "test_get_subsection_entry_list.h"
#include "test_get_entry_name_ending_list.h"
#include "test_load_BC_list.h"


#include "../../include/SolverParameters.h"
#include "../../src/SolverParameters/SolverParameters.cc"

#include "test_LinearSolverParameters.h"
#include "test_NonlinearSolverParameters.h"

#include "../../src/EquationDependencyParser/EquationDependencyParser.cc"
#include "test_EquationDependencyParser.h"
