// List of variables and residual equations for the grain growth example application

// =================================================================================
// Define the variables in the model
// =================================================================================
// The number of variables
#define num_var 10

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define variable_name {"n1","n2","n3","n4","n5","n6","n7","n8","n9","n10"}
#define variable_type {"SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR"}
#define variable_eq_type {"PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define need_val {true,true,true,true,true,true,true,true,true,true}
#define need_grad {true,true,true,true,true,true,true,true,true,true}

// Flags for whether the residual equation has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define need_hess  {false,false,false,false,false,false,false,false,false,false}
#define need_val_residual {true,true,true,true,true,true,true,true,true,true}
#define need_grad_residual {true,true,true,true,true,true,true,true,true,true}

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".

// Mobility
#define MnV 1.0

// Gradient energy coefficient
#define KnV 0.6

// Grain interaction coefficient
#define alpha 1.5

// All residuals defined below

// =================================================================================
// residualRHS
// =================================================================================
// This function calculates the residual equations for each variable. It takes
// "modelVariablesList" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of
// that quadrature point is given by "q_point_loc". The function outputs
// "modelResidualsList", a list of the value and gradient terms of the residual for
// each residual equation. The index for each variable in these lists corresponds to
// the order it is defined at the top of this file (starting at 0).
template <int dim>
void generalizedProblem<dim>::residualRHS(const std::vector<modelVariable<dim>> & modelVariablesList,
												std::vector<modelResidual<dim>> & modelResidualsList,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

dealii::VectorizedArray<double> fnV = constV(0.0);
scalargradType nx;

for (unsigned int i=0; i<num_var; i++){
	fnV = - modelVariablesList[i].scalarValue + modelVariablesList[i].scalarValue*modelVariablesList[i].scalarValue*modelVariablesList[i].scalarValue;
	nx = modelVariablesList[i].scalarGrad;
	for (unsigned int j=0; j<num_var; j++){
		if (i != j){
			fnV += constV(2.0*alpha) * modelVariablesList[i].scalarValue * modelVariablesList[j].scalarValue*modelVariablesList[j].scalarValue;
		}
	}
	modelResidualsList[i].scalarValueResidual = modelVariablesList[i].scalarValue-constV(timeStep*MnV)*fnV;
	modelResidualsList[i].scalarGradResidual = constV(-timeStep*KnV*MnV)*nx;
}


}

// =================================================================================
// residualLHS (needed only if at least one equation is elliptic)
// =================================================================================
// This function calculates the residual equations for the iterative solver for
// elliptic equations.for each variable. It takes "modelVariablesList" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is given
// by "q_point_loc". The function outputs "modelRes", the value and gradient terms of
// for the left-hand-side of the residual equation for the iterative solver. The
// index for each variable in these lists corresponds to the order it is defined at
// the top of this file (starting at 0), not counting variables that have
// "need_val_LHS", "need_grad_LHS", and "need_hess_LHS" all set to "false". If there
// are multiple elliptic equations, conditional statements should be used to ensure
// that the correct residual is being submitted. The index of the field being solved
// can be accessed by "this->currentFieldIndex".
template <int dim>
void generalizedProblem<dim>::residualLHS(const std::vector<modelVariable<dim>> & modelVarList,
		modelResidual<dim> & modelRes,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}

// =================================================================================
// energyDensity (needed only if calcEnergy == true)
// =================================================================================
// This function integrates the free energy density across the computational domain.
// It takes "modelVariablesList" as an input, which is a list of the value and
// derivatives of each of the variables at a specific quadrature point. It also
// takes the mapped quadrature weight, "JxW_value", as an input. The (x,y,z) location
// of the quadrature point is given by "q_point_loc". The weighted value of the
// energy density is added to "energy" variable and the components of the energy
// density are added to the "energy_components" variable (index 0: chemical energy,
// index 1: gradient energy, index 2: elastic energy).
template <int dim>
void generalizedProblem<dim>::energyDensity(const std::vector<modelVariable<dim>> & modelVarList,
											const dealii::VectorizedArray<double> & JxW_value,
											dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {
	scalarvalueType total_energy_density = constV(0.0);


}




