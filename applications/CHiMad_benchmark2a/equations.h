// List of variables and residual equations for the coupled Allen-Cahn/Cahn-Hilliard example application

// =================================================================================
// Define the variables in the model
// =================================================================================
// The number of variables
#define num_var 6

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define variable_name {"c", "mu", "n1", "n2", "n3", "n4"}
#define variable_type {"SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR"}
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
#define McV constV(5.0)

// Allen-Cahn mobility
#define MnV constV(5.0)

// Allen-Cahn gradient energy coefficient
#define KnV constV(3.0)

// Cahn-Hilliard gradient energy coefficient
#define KcV constV(3.0)

//Common coefficient for the double well and structural order parameter interaction terms
#define wV constV(1.0)

//Coefficient order parameter for the structural interaction term only
#define alpha constV(5.0)

// Free energy for each phase and their first and second derivatives
#define faV (constV(2.0)*(c-constV(0.3))*(c-constV(0.3)))
#define facV (constV(4.0)*(c-constV(0.3)))
#define faccV (constV(4.0))
#define fbV (constV(2.0)*(c-constV(0.7))*(c-constV(0.7)))
#define fbcV (constV(4.0)*(c-constV(0.7)))
#define fbccV (constV(4.0))

// Interpolation function and its derivatives
#define hV ( n1*n1*n1*(constV(6.0)*n1*n1-constV(15.0)*n1+constV(10.0)) + n2*n2*n2*(constV(6.0)*n2*n2-constV(15.0)*n2+constV(10.0)) + n3*n3*n3*(constV(6.0)*n3*n3-constV(15.0)*n3+constV(10.0)) + n4*n4*n4*(constV(6.0)*n4*n4-constV(15.0)*n4+constV(10.0)) )
#define hn1V ( n1*n1*(constV(30.0)*n1*n1-constV(60.0)*n1+constV(30.0)) )
#define hn2V ( n2*n2*(constV(30.0)*n2*n2-constV(60.0)*n2+constV(30.0)) )
#define hn3V ( n3*n3*(constV(30.0)*n3*n3-constV(60.0)*n3+constV(30.0)) )
#define hn4V ( n4*n4*(constV(30.0)*n4*n4-constV(60.0)*n4+constV(30.0)) )

//Combined double-well and interaction functions (function g) and its derivatives
//Double-well part
#define gdwV ( n1*n1*(constV(1.0)-n1)*(constV(1.0)-n1) + n2*n2*(constV(1.0)-n2)*(constV(1.0)-n2) +n3*n3*(constV(1.0)-n3)*(constV(1.0)-n3) + n4*n4*(constV(1.0)-n4)*(constV(1.0)-n4) )
//Interaction part
#define gintV ( alpha*(n1*n1*n2*n2 + n1*n1*n3*n3+ n1*n1*n4*n4 + n2*n2*n3*n3 + n2*n2*n4*n4 + n3*n3*n4*n4) )
//Combined function (g)
#define gV ( gdwV + gintV )
//Derivatives
#define dgn1V ( constV(2.0)*n1*(constV(1.0)-n1)*(constV(1.0)-constV(2.0)*n1) + constV(2.0)*alpha*n1*(n2*n2+n3*n3+n4*n4) )
#define dgn2V ( constV(2.0)*n2*(constV(1.0)-n2)*(constV(1.0)-constV(2.0)*n2) + constV(2.0)*alpha*n2*(n1*n1+n3*n3+n4*n4) )
#define dgn3V ( constV(2.0)*n3*(constV(1.0)-n3)*(constV(1.0)-constV(2.0)*n3) + constV(2.0)*alpha*n3*(n1*n1+n2*n2+n4*n4) )
#define dgn4V ( constV(2.0)*n4*(constV(1.0)-n4)*(constV(1.0)-constV(2.0)*n4) + constV(2.0)*alpha*n4*(n1*n1+n2*n2+n3*n3) )

//CHANGES FOR CHAC BENCHMARK PROBLEM UP TO THIS POINT

// Residual equations
#define rn1V ( n1 - constV(timeStep)*MnV*((fbV-faV)*hn1V + wV*dgn1V) )
#define rn2V ( n2 - constV(timeStep)*MnV*((fbV-faV)*hn2V + wV*dgn2V) )
#define rn3V ( n3 - constV(timeStep)*MnV*((fbV-faV)*hn3V + wV*dgn3V) )
#define rn4V ( n4 - constV(timeStep)*MnV*((fbV-faV)*hn4V + wV*dgn4V) )
#define rn1xV ( constV(-timeStep)*KnV*MnV*n1x )
#define rn2xV ( constV(-timeStep)*KnV*MnV*n2x )
#define rn3xV ( constV(-timeStep)*KnV*MnV*n3x )
#define rn4xV ( constV(-timeStep)*KnV*MnV*n4x )
#define rmuV ( (constV(1.0)-hV)*facV+hV*fbcV )
#define rmuxV ( KcV*cx )
#define rcV ( c )
#define rcxV ( constV(-timeStep)*McV*mux )


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

// The chemical potential and its derivatives (names here should match those in the macros above)
scalargradType mux = modelVariablesList[1].scalarGrad;

// The order parameters and their derivatives (names here should match those in the macros above)
scalarvalueType n1 = modelVariablesList[2].scalarValue;
scalargradType n1x = modelVariablesList[2].scalarGrad;
scalarvalueType n2 = modelVariablesList[3].scalarValue;
scalargradType n2x = modelVariablesList[3].scalarGrad;
scalarvalueType n3 = modelVariablesList[4].scalarValue;
scalargradType n3x = modelVariablesList[4].scalarGrad;
scalarvalueType n4 = modelVariablesList[5].scalarValue;
scalargradType n4x = modelVariablesList[5].scalarGrad;

// Residuals for the equation to evolve the concentration (names here should match those in the macros above)
modelResidualsList[0].scalarValueResidual = rcV;
modelResidualsList[0].scalarGradResidual = rcxV;
modelResidualsList[1].scalarValueResidual = rmuV;
modelResidualsList[1].scalarGradResidual = rmuxV;

// Residuals for the equation to evolve the order parameters (names here should match those in the macros above)
modelResidualsList[2].scalarValueResidual = rn1V;
modelResidualsList[2].scalarGradResidual = rn1xV;
modelResidualsList[3].scalarValueResidual = rn2V;
modelResidualsList[3].scalarGradResidual = rn2xV;
modelResidualsList[4].scalarValueResidual = rn3V;
modelResidualsList[4].scalarGradResidual = rn3xV;
modelResidualsList[5].scalarValueResidual = rn4V;
modelResidualsList[5].scalarGradResidual = rn4xV;

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
scalarvalueType n1 = modelVariablesList[2].scalarValue;
scalargradType n1x = modelVariablesList[2].scalarGrad;
scalarvalueType n2 = modelVariablesList[3].scalarValue;
scalargradType n2x = modelVariablesList[3].scalarGrad;
scalarvalueType n3 = modelVariablesList[4].scalarValue;
scalargradType n3x = modelVariablesList[4].scalarGrad;
scalarvalueType n4 = modelVariablesList[5].scalarValue;
scalargradType n4x = modelVariablesList[5].scalarGrad;

// The homogenous free energy
scalarvalueType f_chem = (constV(1.0)-hV)*faV + hV*fbV + wV*gV;

// The gradient free energy
scalarvalueType f_grad = constV(0.5)*KnV*(n1x*n1x+n2x*n2x+n3x*n3x+n4x*n4x) + constV(0.5)*KcV*cx*cx;

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




