// List of variables and residual equations for the coupled Allen-Cahn/Cahn-Hilliard example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
<<<<<<< HEAD
// The number of variables
#define num_var 6

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define variable_name {"c", "mu", "n1", "n2", "n3", "n4"}
#define variable_type {"SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR"}
#define variable_eq_type {"PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define need_val {true, true, true, true, true, true}
#define need_grad {true, true, true, true, true, true}
#define need_hess {false, false, false, false, false, false}

// Flags for whether the residual equation has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define need_val_residual {true, true, true, true, true, true}
#define need_grad_residual {true, true, true, true, true, true}
=======
void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"c");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,PARABOLIC);

	set_need_value					(0,true);
	set_need_gradient				(0,true);
	set_need_hessian				(0,false);

	set_need_value_residual_term	(0,true);
	set_need_gradient_residual_term	(0,true);

    // Variable 1
	set_variable_name				(1,"n");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,PARABOLIC);

	set_need_value					(1,true);
	set_need_gradient				(1,true);
	set_need_hessian				(1,false);

	set_need_value_residual_term	(1,true);
	set_need_gradient_residual_term	(1,true);
}
>>>>>>> 0147129c14b77ac7320a8032fc57a37558811c95

// =================================================================================
// residualRHS
// =================================================================================
// This function calculates the residual equations for each variable. It takes
// "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of
// that quadrature point is given by "q_point_loc". The function outputs residuals
// to variable_list. The index for each variable in this list corresponds to
// the index given at the top of this file.

<<<<<<< HEAD
// Cahn-Hilliard mobility
#define McV 5.0

// Allen-Cahn mobility
#define MnV 5.0

// Allen-Cahn gradient energy coefficient
#define KnV 3.0

// Cahn-Hilliard gradient energy coefficient
#define KcV 3.0

//Common coefficient for the double well and structural order parameter interaction terms
#define wV 1.0

//Coefficient order parameter for the structural interaction term only
#define alpha 5.0

// Free energy for each phase and their first and second derivatives
#define faV (2.0*(c-0.3)*(c-0.3))
#define facV (4.0*(c-0.3))
#define faccV (4.0)
#define fbV (2.0*(c-0.7)*(c-0.7))
#define fbcV (4.0*(c-0.7))
#define fbccV (4.0)

// Interpolation function and its derivatives
#define hV ( n1*n1*n1*(6.0*n1*n1-15.0*n1+10.0) + n2*n2*n2*(6.0*n2*n2-15.0*n2+10.0) + n3*n3*n3*(6.0*n3*n3-15.0*n3+10.0) + n4*n4*n4*(6.0*n4*n4-15.0*n4+10.0) )
#define hn1V ( n1*n1*(30.0*n1*n1-60.0*n1+30.0) )
#define hn2V ( n2*n2*(30.0*n2*n2-60.0*n2+30.0) )
#define hn3V ( n3*n3*(30.0*n3*n3-60.0*n3+30.0) )
#define hn4V ( n4*n4*(30.0*n4*n4-60.0*n4+30.0) )

//Combined double-well and interaction functions (function g) and its derivatives
//Double-well part
#define gdwV ( n1*n1*(1.0-n1)*(1.0-n1) + n2*n2*(1.0-n2)*(1.0-n2) +n3*n3*(1.0-n3)*(1.0-n3) + n4*n4*(1.0-n4)*(1.0-n4) )
//Interaction part
#define gintV ( 2.0*alpha*(n1*n1*n2*n2 + n1*n1*n3*n3+ n1*n1*n4*n4 + n2*n2*n3*n3 + n2*n2*n4*n4 + n3*n3*n4*n4) )
//Combined function (g)
#define gV ( gdwV + gintV )
//Derivatives
#define dgn1V ( 2.0*n1*(1.0-n1)*(1.0-2.0*n1) + 4.0*alpha*n1*(n2*n2+n3*n3+n4*n4) )
#define dgn2V ( 2.0*n2*(1.0-n2)*(1.0-2.0*n2) + 4.0*alpha*n2*(n1*n1+n3*n3+n4*n4) )
#define dgn3V ( 2.0*n3*(1.0-n3)*(1.0-2.0*n3) + 4.0*alpha*n3*(n1*n1+n2*n2+n4*n4) )
#define dgn4V ( 2.0*n4*(1.0-n4)*(1.0-2.0*n4) + 4.0*alpha*n4*(n1*n1+n2*n2+n3*n3) )

//CHANGES FOR CHAC BENCHMARK PROBLEM UP TO THIS POINT

// Residual equations
#define rn1V ( n1 - constV(timeStep*MnV)*((fbV-faV)*hn1V+constV(wV)*dgn1V) )
#define rn2V ( n2 - constV(timeStep*MnV)*((fbV-faV)*hn2V+constV(wV)*dgn2V) )
#define rn3V ( n3 - constV(timeStep*MnV)*((fbV-faV)*hn3V+constV(wV)*dgn3V) )
#define rn4V ( n4 - constV(timeStep*MnV)*((fbV-faV)*hn4V+constV(wV)*dgn4V) )
#define rn1xV ( constV(-timeStep*KnV*MnV)*n1x )
#define rn2xV ( constV(-timeStep*KnV*MnV)*n2x )
#define rn3xV ( constV(-timeStep*KnV*MnV)*n3x )
#define rn4xV ( constV(-timeStep*KnV*MnV)*n4x )
#define rmuV ( (constV(1.0)-hV)*facV+hV*fbcV )
#define rmuxV ( constV(KcV)*cx )
#define rcV ( c )
#define rcxV ( constV(-timeStep*McV)*mux )


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
scalarvalueType mu = modelVariablesList[1].scalarValue;
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
=======
template <int dim, int degree>
void customPDE<dim,degree>::residualRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//c
scalarvalueType c = variable_list.get_scalar_value(0);
scalargradType cx = variable_list.get_scalar_gradient(0);

//n
scalarvalueType n = variable_list.get_scalar_value(1);
scalargradType nx = variable_list.get_scalar_gradient(1);

// Free energy for each phase and their first and second derivatives
scalarvalueType faV = (-1.6704-4.776*c+5.1622*c*c-2.7375*c*c*c+1.3687*c*c*c*c);
scalarvalueType facV = (-4.776 + 10.3244*c - 8.2125*c*c + 5.4748*c*c*c);
scalarvalueType faccV = (10.3244-16.425*c+16.4244*c*c);
scalarvalueType fbV = (5.0*c*c-5.9746*c-1.5924);
scalarvalueType fbcV = (10.0*c-5.9746);
scalarvalueType fbccV = constV(10.0);

// Interpolation function and its derivative
scalarvalueType hV = (10.0*n*n*n-15.0*n*n*n*n+6.0*n*n*n*n*n);
scalarvalueType hnV = (30.0*n*n-60.0*n*n*n+30.0*n*n*n*n);

// Residual equations
scalargradType muxV = ( cx*((1.0-hV)*faccV+hV*fbccV) + nx*((fbcV-facV)*hnV) );
scalarvalueType rcV = c;
scalargradType rcxV = (constV(-McV*userInputs.dtValue)*muxV);
scalarvalueType rnV = (n-constV(userInputs.dtValue*MnV)*(fbV-faV)*hnV);
scalargradType rnxV = (constV(-userInputs.dtValue*KnV*MnV)*nx);

// Residuals for the equation to evolve the concentration (names here should match those in the macros above)
variable_list.set_scalar_value_residual_term(0,rcV);
variable_list.set_scalar_gradient_residual_term(0,rcxV);

// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
variable_list.set_scalar_value_residual_term(1,rnV);
variable_list.set_scalar_gradient_residual_term(1,rnxV);
>>>>>>> 0147129c14b77ac7320a8032fc57a37558811c95

}

// =================================================================================
// residualLHS (needed only if at least one equation is elliptic)
// =================================================================================
// This function calculates the residual equations for the iterative solver for
// elliptic equations.for each variable. It takes "variable_list" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is given
// by "q_point_loc". The function outputs residual terms to "variable_list"
// for the left-hand-side of the residual equation for the iterative solver. The
// index for each variable in this list corresponds to
// the index given at the top of this file. If there are multiple elliptic equations,
// conditional statements should be used to ensure that the correct residual is
// being submitted. The index of the field being solved can be accessed by
// "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::residualLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
<<<<<<< HEAD

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
scalarvalueType f_chem = (constV(1.0)-hV)*faV + hV*fbV + constV(wV)*gV;

// The gradient free energy
scalarvalueType f_grad = constV(0.5*KnV)*(n1x*n1x+n2x*n2x+n3x*n3x+n4x*n4x) + constV(0.5*KcV)*cx*cx;

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
=======
>>>>>>> 0147129c14b77ac7320a8032fc57a37558811c95
}
