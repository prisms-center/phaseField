// List of variables and residual equations for the coupled Allen-Cahn/Cahn-Hilliard example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"c");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_residual_term_RHS(0, "c");
    set_dependencies_gradient_residual_term_RHS(0, "n,grad(c)");

    // Variable 1
	set_variable_name				(1,"n");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_residual_term_RHS(1, "c,n");
    set_dependencies_gradient_residual_term_RHS(1, "grad(n)");
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
void customPDE<dim,degree>::residualExplicitRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//c
scalarvalueType c = variable_list.get_scalar_value(0);
scalargradType cx = variable_list.get_scalar_gradient(0);

//n
scalarvalueType n = variable_list.get_scalar_value(1);
scalargradType nx = variable_list.get_scalar_gradient(1);

// Free energy for each phase and their first and second derivatives
scalarvalueType fa = constV(2.0)*c*c;
scalarvalueType fac = constV(4.0)*c;
scalarvalueType facc = constV(4.0);
scalarvalueType fb = constV(2.0)*(c*c - 2.0*c + constV(1.0));
scalarvalueType fbc = constV(4.0)*(c - 1.0);
scalarvalueType fbcc = constV(4.0);

// Interpolation function and its derivative
scalarvalueType h = (3.0*n*n-2.0*n*n*n);
scalarvalueType hn = (6.0*n-6.0*n*n);

// Residual equations
scalargradType mux = ( cx*((1.0-h)*facc+h*fbcc) + nx*((fbc-fac)*hn) );
scalarvalueType rc = c;
scalargradType rcx = (constV(-Mc*userInputs.dtValue)*mux);
scalarvalueType rn = (n-constV(userInputs.dtValue*Mn)*(fb-fa)*hn);
scalargradType rnx = (constV(-userInputs.dtValue*Kn*Mn)*nx);

// Residuals for the equation to evolve the concentration
variable_list.set_scalar_value_residual_term(0,rc);
variable_list.set_scalar_gradient_residual_term(0,rcx);

// Residuals for the equation to evolve the order parameter
variable_list.set_scalar_value_residual_term(1,rn);
variable_list.set_scalar_gradient_residual_term(1,rnx);

}

template <int dim, int degree>
void customPDE<dim,degree>::residualNonexplicitRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

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
