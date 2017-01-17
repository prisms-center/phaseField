// List of variables and residual equations for the coupled Allen-Cahn/Cahn-Hilliard example application

// =================================================================================
// Define the variables in the model
// =================================================================================
// The number of variables
#define num_var 5

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define variable_name {"c", "mu","n1", "n2", "n3", "n4"}
#define variable_type {"SCALAR","SCALAR", "SCALAR","SCALAR","SCALAR","SCALAR"}
#define variable_eq_type {"PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define need_val {true, false, true, true, true, true}
#define need_grad {true, true, true, true, true, true}
#define need_hess {false, false, false, false, false, false}

// Flags for whether the residual equation has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define need_val_residual {true, true, true, true, true, true}
#define need_grad_residual {true, true, true, true, true, true}

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".

// Cahn-Hilliard mobility
#define McV 5.0

// Allen-Cahn mobility
#define MnV 5.0

// Allen-Cahn gradient energy coefficient
#define KnV 3.0

// Cahn-Hilliard gradient energy coefficient
#define KcV 3.0

//Common coefficient for the double well and structural order parameter interaction terms
#define w 1.0

//Coefficient order parameter for the structural interaction term only
#define alpha 5.0

// Free energy for each phase and their first and second derivatives
#define faV (2.0*(c-0.3)*(c-0.3))
#define facV (4.0*(c-0.3))
#define faccV (4.0)
#define fbV (2.0*(c-0.7)*(c-0.7))
#define fbcV (4.0*(c-0.3))
#define fbccV (4.0)

// Interpolation function and its derivative
#define hV ( n1*n1*n1*(6.0*n1*n1-15.0*n1+10.0) + n2*n2*n2*(6.0*n2*n2-15.0*n2+10.0) + n3*n3*n3*(6.0*n3*n3-15.0*n3+10.0) + n4*n4*n4*(6.0*n4*n4-15.0*n4+10.0) )
#define hn1V ( n1*n1*(30.0*n1*n1-60*n1+30.0) )
#define hn2V ( n2*n2*(30.0*n2*n2-60*n2+30.0) )
#define hn3V ( n3*n3*(30.0*n3*n3-60*n3+30.0) )
#define hn4V ( n4*n4*(30.0*n4*n4-60*n4+30.0) )

//Derivative of combined (function g) double-well and interaction functions

#define gn1 ( 2.0*n1*(1.0-n1)*(1.0-2.0*n1) + 2.0*alpha*n1*(n2*n2+n3*n3+n4*n4) )
#define gn2 ( 2.0*n2*(1.0-n2)*(1.0-2.0*n2) + 2.0*alpha*n2*(n1*n1+n3*n3+n4*n4) )
#define gn3 ( 2.0*n3*(1.0-n3)*(1.0-2.0*n3) + 2.0*alpha*n3*(n1*n1+n2*n2+n4*n4) )
#define gn4 ( 2.0*n4*(1.0-n4)*(1.0-2.0*n4) + 2.0*alpha*n4*(n1*n1+n2*n2+n3*n3) )

//CHANGES FOR CHAC BENCHMARK PROBLEM UP TO THIS POINT

// Residual equations
#define rmuV ( hV * (fbcV - facV) )
#define rmuxV ( constV(KcV) * cx  )
#define rcV   (c)
#define rcxV  (constV(-McV*timeStep)*mux)
#define rn1V  (n1-constV(timeStep*MnV)*(fbV-faV)*hnV)
#define rn1xV (constV(-timeStep*KnV*MnV)*nx)
#define rnV  (n-constV(timeStep*MnV)*(fbV-faV)*hnV)
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
template <int dim>
void generalizedProblem<dim>::residualRHS(const std::vector<modelVariable<dim> > & modelVariablesList,
												std::vector<modelResidual<dim> > & modelResidualsList,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// The concentration and its derivatives (names here should match those in the macros above)
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

// The order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n = modelVariablesList[1].scalarValue;
scalargradType nx = modelVariablesList[1].scalarGrad;

// Residuals for the equation to evolve the concentration (names here should match those in the macros above)
modelResidualsList[0].scalarValueResidual = rcV;
modelResidualsList[0].scalarGradResidual = rcxV;

// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
modelResidualsList[1].scalarValueResidual = rnV;
modelResidualsList[1].scalarGradResidual = rnxV;

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
void generalizedProblem<dim>::residualLHS(const std::vector<modelVariable<dim> > & modelVariablesList,
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
void generalizedProblem<dim>::energyDensity(const std::vector<modelVariable<dim> > & modelVariablesList,
											const dealii::VectorizedArray<double> & JxW_value,
											dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {

// The concentration and its derivatives (names here should match those in the macros above)
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

// The order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n = modelVariablesList[1].scalarValue;
scalargradType nx = modelVariablesList[1].scalarGrad;

// The homogenous free energy
scalarvalueType f_chem = (constV(1.0)-hV)*faV + hV*fbV;

// The gradient free energy
scalarvalueType f_grad = constV(0.5*KnV)*nx*nx;

// The total free energy
scalarvalueType total_energy_density;
total_energy_density = f_chem + f_grad;

// Loop to step through each element of the vectorized arrays. Working with deal.ii
// developers to see if there is a more elegant way to do this.
assembler_lock.acquire ();
for (unsigned i=0; i<c.n_array_elements;i++){
  if (c[i] > 1.0e-10){
	  this->energy+=total_energy_density[i]*JxW_value[i];
	  this->energy_components[0]+= f_chem[i]*JxW_value[i];
	  this->energy_components[1]+= f_grad[i]*JxW_value[i];
  }
}
assembler_lock.release ();
}




