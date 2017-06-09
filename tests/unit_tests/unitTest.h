#include "../../include/dealIIheaders.h"
#include <iostream>

#define problemDIM 2
#define finiteElementDegree 1
#define vectorgradType dealii::Tensor<2, dim, dealii::VectorizedArray<double> >
#define typeScalar dealii::FEEvaluation<dim,finiteElementDegree,finiteElementDegree+1,1,double>
#define typeVector dealii::FEEvaluation<dim,finiteElementDegree,finiteElementDegree+1,dim,double>

//define test variables for the tests
#define subdivisionsX 10
#define subdivisionsY 10
#define subdivisionsZ 10

#define numOutputs 10
#define timeStep 1.0e-3
#define timeFinal 20.0
#define timeIncrements 20000

// =================================================================================
// Define the variables in the model
// =================================================================================
// The number of variables
#define num_var 1

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define variable_name {"n"}
#define variable_type {"SCALAR"}
#define variable_eq_type {"PARABOLIC"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define need_val {true}
#define need_grad {true}
#define need_hess  {false}

// Flags for whether the residual equation has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define need_val_residual {true}
#define need_grad_residual {true}



//define data type
template <int dim>
void computeStress(const dealii::Table<2, double>& CIJ, const dealii::VectorizedArray<double> ux[][dim], const dealii::VectorizedArray<double> R[][dim]);

#include "../../include/matrixFreePDE.h"
#include "../../include/postprocessor.h"

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
#include "../../src/matrixfree/calcFreeEnergy.cc"
#include "../../src/matrixfree/computeIntegral.cc"
#include "../../src/matrixfree/getOutputTimeSteps.cc"

#include "../../src/inputFileReader/inputFileReader.cc"

#include "../../src/models/mechanics/computeStress.h"

#define POSTPROCESSOR_TEMPLATES
#include "postprocess.h"
#include "../../src/postprocessor/postprocessor.cc"

template <int dim, typename T>
class unitTest
{
	public:
	bool test_computeInvM(int argc, char **argv, userInputParameters<dim>);
	bool test_outputResults(int argc, char **argv, userInputParameters<dim> userInputs);
	bool test_computeStress();
	void assignCIJSize(dealii::VectorizedArray<double> CIJ[2*dim-1+dim/3][2*dim-1+dim/3]);
	void assignCIJSize(dealii::Table<2, double> &CIJ);
	bool test_getOutputTimeSteps(std::string,unsigned int, std::vector<unsigned int>, userInputParameters<dim> userInputs);
	bool test_setRigidBodyModeConstraints(std::vector<int>, userInputParameters<dim> userInputs);
	bool test_parse_line();
	bool test_get_subsection_entry_list();
	bool test_load_BC_list();
};


#include "test_invM.h"
#include "test_outputResults.h"
#include "test_computeStress.h"
#include "test_getOutputTimeSteps.h"
#include "test_setRigidBodyModeConstraints.h"
#include "test_parse_line.h"
#include "test_get_subsection_entry_list.h"
#include "test_load_BC_list.h"
