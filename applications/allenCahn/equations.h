// List of variables and residual equations for the coupled Allen-Cahn example application

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
#define KnV 2.0

// Free energy and its derivative
#define fV (n*n*n*n - 2.0*n*n*n + n*n)
#define fnV (4.0*n*(n-1.0)*(n-0.5))

// Residual equations
#define rnV  (n-constV(timeStep*MnV)*fnV)
#define rnxV (constV(-timeStep*KnV*MnV)*nx)

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
template <int dim, int degree>
void customPDE<dim,degree>::residualRHS(const std::vector<modelVariable<dim> > & modelVariablesList,
												std::vector<modelResidual<dim> > & modelResidualsList,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// The order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n = modelVariablesList[0].scalarValue;
scalargradType nx = modelVariablesList[0].scalarGrad;


// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
modelResidualsList[0].scalarValueResidual = rnV;
modelResidualsList[0].scalarGradResidual = rnxV;

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
template <int dim, int degree>
void customPDE<dim,degree>::residualLHS(const std::vector<modelVariable<dim> > & modelVarList,
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
template <int dim, int degree>
void customPDE<dim,degree>::energyDensity(const std::vector<modelVariable<dim> > & modelVarList,
											const dealii::VectorizedArray<double> & JxW_value,
											dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {

scalarvalueType total_energy_density = constV(0.0);

// The order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n = modelVarList[0].scalarValue;
scalargradType nx = modelVarList[0].scalarGrad;

// The homogenous free energy
scalarvalueType f_chem = fV;

// The gradient free energy
scalarvalueType f_grad = constV(0.0);

for (int i=0; i<dim; i++){
  for (int j=0; j<dim; j++){
	  f_grad += constV(0.5*KnV)*nx[i]*nx[j];
  }
}

// The total free energy
total_energy_density = f_chem + f_grad;

// Loop to step through each element of the vectorized arrays. Working with deal.ii
// developers to see if there is a more elegant way to do this.
this->assembler_lock.acquire ();
for (unsigned i=0; i<n.n_array_elements;i++){
  // For some reason, some of the values in this loop
  if (n[i] > 1.0e-10){
	  this->energy+=total_energy_density[i]*JxW_value[i];
	  this->energy_components[0]+= f_chem[i]*JxW_value[i];
	  this->energy_components[1]+= f_grad[i]*JxW_value[i];
  }
}
this->assembler_lock.release ();
}
