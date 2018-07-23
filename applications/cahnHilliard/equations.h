// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"c");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_residual_term_RHS(0, "c");
    set_dependencies_gradient_residual_term_RHS(0, "grad(mu)");

	// Variable 1
	set_variable_name				(1,"mu");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,AUXILIARY);

    set_dependencies_value_residual_term_RHS(1, "c");
    set_dependencies_gradient_residual_term_RHS(1, "grad(c)");

}

// =================================================================================
// explicitEquationRHS
// =================================================================================
// This function calculates the residual equations for each variable. It takes
// "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of
// that quadrature point is given by "q_point_loc". The function outputs residuals
// to variable_list. The index for each variable in this list corresponds to
// the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// The concentration and its derivatives
scalarvalueType c = variable_list.get_scalar_value(0);

// The chemical potential and its derivatives
scalargradType mux = variable_list.get_scalar_gradient(1);

// Parameters in the residual equations and expressions for the residual equations
// can be set here.


// The residuals
scalarvalueType rcV = c;
scalargradType rcxV = constV(-McV*userInputs.dtValue)*mux;

// Residuals for the equation to evolve the concentration
variable_list.set_scalar_value_residual_term(0,rcV);
variable_list.set_scalar_gradient_residual_term(0,rcxV);

}

// =================================================================================
// nonExplicitEquationRHS
// =================================================================================
// This function calculates the residual equations for each variable. It takes
// "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of
// that quadrature point is given by "q_point_loc". The function outputs residuals
// to variable_list. The index for each variable in this list corresponds to
// the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

 // The concentration and its derivatives
 scalarvalueType c = variable_list.get_scalar_value(0);
 scalargradType cx = variable_list.get_scalar_gradient(0);

 // Parameters in the residual equations and expressions for the residual equations
 // can be set here.

 // The derivative of the local free energy
 scalarvalueType fcV = 4.0*c*(c-1.0)*(c-0.5);

 // The residuals
 scalarvalueType rmuV = fcV;
 scalargradType rmuxV = constV(KcV)*cx;

 // Residuals for the equation to evolve the chemical potential
 variable_list.set_scalar_value_residual_term(1,rmuV);
 variable_list.set_scalar_gradient_residual_term(1,rmuxV);


}

// =================================================================================
// equationLHS (needed only if at least one equation is elliptic)
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
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}
