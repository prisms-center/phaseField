// List of residual equations for the coupled Allen-Cahn example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"n");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,PARABOLIC);

	set_need_value					(0,true);
	set_need_gradient				(0,true);
	set_need_hessian					(0,false);

	set_need_value_residual_term		(0,true);
	set_need_gradient_residual_term	(0,true);
}

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

// The order parameter and its derivatives 
scalarvalueType n = variable_list.get_scalar_value(0); //modelVariablesList[0].scalarValue;
scalargradType nx = variable_list.get_scalar_gradient(0); //modelVariablesList[0].scalarGrad;

// Parameters in the residual equations and expressions for the residual equations
// can be set here.
scalarvalueType fnV = (4.0*n*(n-1.0)*(n-0.5));
scalarvalueType rnV = (n-constV(userInputs.dtValue*MnV)*fnV);
scalargradType rnxV = (-constV(userInputs.dtValue*KnV*MnV)*nx);

// Residuals for the equation to evolve the order parameter 
variable_list.set_scalar_value_residual_term(0,rnV);
variable_list.set_scalar_gradient_residual_term(0,rnxV);

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
