// List of variables and residual equations for the coupled Allen-Cahn/Cahn-Hilliard example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"T");
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
	set_need_gradient_residual_term	(1,false);

	// Variable 2
	set_variable_name				(2,"mu");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,PARABOLIC);

	set_need_value					(2,true);
	set_need_gradient				(2,false);
	set_need_hessian				(2,false);

	set_need_value_residual_term	(2,true);
	set_need_gradient_residual_term	(2,true);
}

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".


// Free energy expression (or rather, its derivative)



// =================================================================================
// residualRHS
// =================================================================================
// This function calculates the residual equations for each variable. It takes
// "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of
// that quadrature point is given by "q_point_loc". The function outputs residuals
// to variable_list. The index for each variable in this list corresponds to
// the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::residualRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// The temperature and its derivatives (names here should match those in the macros above)
scalarvalueType T = variable_list.get_scalar_value(0);
scalargradType Tx = variable_list.get_scalar_gradient(0);

// The order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n = variable_list.get_scalar_value(1);
scalargradType nx = variable_list.get_scalar_gradient(1);

// The order parameter chemical potential and its derivatives (names here should match those in the macros above)
scalarvalueType mu = variable_list.get_scalar_value(2);

double lambdaV = (DV/0.6267/W0/W0);

// Derivative of the free energy density with respect to n
scalarvalueType fnV = (-(n-constV(lambdaV)*T*(constV(1.0)-n*n))*(constV(1.0)-n*n));

scalarvalueType theta;

for (unsigned i=0; i< n.n_array_elements;i++){
	theta[i] = std::atan2(nx[1][i],nx[0][i]);
}

// Anisotropic gradient energy coefficient, its derivative and square
scalarvalueType WV = (constV(W0)*(constV(1.0)+constV(epsilonM)*std::cos(constV(mult)*(theta-constV(theta0)))));
scalarvalueType WTV = (constV(W0)*(-constV(epsilonM)*constV(mult)*std::sin(constV(mult)*(theta-constV(theta0)))));
scalarvalueType tauV = (WV*WV);



scalargradType aniso;
aniso[0] = -WV*WV*nx[0]+WV*WTV*nx[1];
aniso[1] = -WV*WV*nx[1]-WV*WTV*nx[0];

// Define required residuals
scalarvalueType rTV = (T-constV(0.5)*mu*constV(userInputs.dtValue)/tauV);
scalargradType rTxV = (constV(-DV*userInputs.dtValue)*Tx);
scalarvalueType rnV = (n-constV(userInputs.dtValue)*mu/tauV);
scalarvalueType rmuV = (fnV);
scalargradType rmuxV = (-aniso);

// Residuals for the equation to evolve the concentration (names here should match those in the macros above)
variable_list.set_scalar_value_residual_term(0,rTV);
variable_list.set_scalar_gradient_residual_term(0,rTxV);

// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
variable_list.set_scalar_value_residual_term(1,rnV);

// Residuals for the equation to evolve the order parameter chemical potential (names here should match those in the macros above)
variable_list.set_scalar_value_residual_term(2,rmuV);
variable_list.set_scalar_gradient_residual_term(2,rmuxV);

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
}
