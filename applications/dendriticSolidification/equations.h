// List of variables and residual equations for the coupled Allen-Cahn/Cahn-Hilliard example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"u");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,PARABOLIC);

	set_need_value					(0,true);
	set_need_gradient				(0,true);
	set_need_hessian				(0,false);

	set_need_value_residual_term	(0,true);
	set_need_gradient_residual_term	(0,true);

    // Variable 1
	set_variable_name				(1,"phi");
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

// The temperature and its derivatives
scalarvalueType u = variable_list.get_scalar_value(0);
scalargradType ux = variable_list.get_scalar_gradient(0);

// The order parameter and its derivatives
scalarvalueType phi = variable_list.get_scalar_value(1);
scalargradType phix = variable_list.get_scalar_gradient(1);

// The order parameter chemical potential and its derivatives
scalarvalueType mu = variable_list.get_scalar_value(2);

// The coupling constant, determined from solvability theory
double lambda = (D/0.6267/W0/W0);

// Derivative of the free energy density with respect to phi
scalarvalueType f_phi = -(phi-constV(lambda)*u*(constV(1.0)-phi*phi))*(constV(1.0)-phi*phi);

// The azimuthal angle
scalarvalueType theta;
for (unsigned i=0; i< phi.n_array_elements;i++){
	theta[i] = std::atan2(phix[1][i],phix[0][i]);
}

// Anisotropic gradient energy coefficient, its derivative and square
scalarvalueType W = constV(W0)*(constV(1.0)+constV(epsilonM)*std::cos(constV(mult)*(theta-constV(theta0))));
scalarvalueType W_theta = constV(-W0)*(constV(epsilonM)*constV(mult)*std::sin(constV(mult)*(theta-constV(theta0))));
scalarvalueType tau = W/constV(W0);

// The anisotropy term that enters in to the residual equation for mu
scalargradType aniso;
aniso[0] = W*W*phix[0]-W*W_theta*phix[1];
aniso[1] = W*W*phix[1]+W*W_theta*phix[0];

// Define required residuals
scalarvalueType ruV = (u+constV(0.5)*mu*constV(userInputs.dtValue)/tau);
scalargradType ruxV = (constV(-D*userInputs.dtValue)*ux);
scalarvalueType rphiV = (phi+constV(userInputs.dtValue)*mu/tau);
scalarvalueType rmuV = (-f_phi);
scalargradType rmuxV = (-aniso);

// Residuals for the equation to evolve the concentration
variable_list.set_scalar_value_residual_term(0,ruV);
variable_list.set_scalar_gradient_residual_term(0,ruxV);

// Residuals for the equation to evolve the order parameter
variable_list.set_scalar_value_residual_term(1,rphiV);

// Residuals for the equation to evolve the order parameter chemical potential
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
