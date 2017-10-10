// List of variables and residual equations for the coupled Allen-Cahn/Cahn-Hilliard example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
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
	set_variable_name				(1,"mu");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,PARABOLIC);

	set_need_value					(1,false);
	set_need_gradient				(1,true);
	set_need_hessian				(1,false);

	set_need_value_residual_term	(1,true);
	set_need_gradient_residual_term	(1,true);

	// Variable 2
	set_variable_name				(2,"n1");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,PARABOLIC);

	set_need_value					(2,true);
	set_need_gradient				(2,true);
	set_need_hessian				(2,false);

	set_need_value_residual_term	(2,true);
	set_need_gradient_residual_term	(2,true);

	// Variable 3
	set_variable_name				(3,"n2");
	set_variable_type				(3,SCALAR);
	set_variable_equation_type		(3,PARABOLIC);

	set_need_value					(3,true);
	set_need_gradient				(3,true);
	set_need_hessian				(3,false);

	set_need_value_residual_term	(3,true);
	set_need_gradient_residual_term	(3,true);

	// Variable 4
	set_variable_name				(4,"n3");
	set_variable_type				(4,SCALAR);
	set_variable_equation_type		(4,PARABOLIC);

	set_need_value					(4,true);
	set_need_gradient				(4,true);
	set_need_hessian				(4,false);

	set_need_value_residual_term	(4,true);
	set_need_gradient_residual_term	(4,true);

	// Variable 5
	set_variable_name				(5,"n4");
	set_variable_type				(5,SCALAR);
	set_variable_equation_type		(5,PARABOLIC);

	set_need_value					(5,true);
	set_need_gradient				(5,true);
	set_need_hessian				(5,false);

	set_need_value_residual_term	(5,true);
	set_need_gradient_residual_term	(5,true);
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

// The concentration and its derivatives 
scalarvalueType c = variable_list.get_scalar_value(0);
scalargradType cx = variable_list.get_scalar_gradient(0);

// The chemical potential and its derivatives 
scalargradType mux = variable_list.get_scalar_gradient(1);

// The order parameters and their derivatives 
scalarvalueType n1 = variable_list.get_scalar_value(2);
scalargradType n1x = variable_list.get_scalar_gradient(2);
scalarvalueType n2 = variable_list.get_scalar_value(3);
scalargradType n2x = variable_list.get_scalar_gradient(3);
scalarvalueType n3 = variable_list.get_scalar_value(4);
scalargradType n3x = variable_list.get_scalar_gradient(4);
scalarvalueType n4 = variable_list.get_scalar_value(5);
scalargradType n4x = variable_list.get_scalar_gradient(5);

// Free energy for each phase and their first and second derivatives
scalarvalueType faV = (constV(2.0)*(c-constV(0.3))*(c-constV(0.3)));
scalarvalueType facV = (constV(4.0)*(c-constV(0.3)));
scalarvalueType faccV = (constV(4.0));
scalarvalueType fbV = (constV(2.0)*(c-constV(0.7))*(c-constV(0.7)));
scalarvalueType fbcV = (constV(4.0)*(c-constV(0.7)));
scalarvalueType fbccV = (constV(4.0));

// Interpolation function and its derivatives
scalarvalueType hV = ( n1*n1*n1*(constV(6.0)*n1*n1-constV(15.0)*n1+constV(10.0)) + n2*n2*n2*(constV(6.0)*n2*n2-constV(15.0)*n2+constV(10.0)) + n3*n3*n3*(constV(6.0)*n3*n3-constV(15.0)*n3+constV(10.0)) + n4*n4*n4*(constV(6.0)*n4*n4-constV(15.0)*n4+constV(10.0)) );
scalarvalueType hn1V = ( n1*n1*(constV(30.0)*n1*n1-constV(60.0)*n1+constV(30.0)) );
scalarvalueType hn2V = ( n2*n2*(constV(30.0)*n2*n2-constV(60.0)*n2+constV(30.0)) );
scalarvalueType hn3V = ( n3*n3*(constV(30.0)*n3*n3-constV(60.0)*n3+constV(30.0)) );
scalarvalueType hn4V = ( n4*n4*(constV(30.0)*n4*n4-constV(60.0)*n4+constV(30.0)) );

//Combined double-well and interaction functions (function g) and its derivatives
//Double-well part
scalarvalueType gdwV = ( n1*n1*(constV(1.0)-n1)*(constV(1.0)-n1) + n2*n2*(constV(1.0)-n2)*(constV(1.0)-n2) +n3*n3*(constV(1.0)-n3)*(constV(1.0)-n3) + n4*n4*(constV(1.0)-n4)*(constV(1.0)-n4) );
//Interaction part
scalarvalueType gintV = ( alpha*(n1*n1*n2*n2 + n1*n1*n3*n3+ n1*n1*n4*n4 + n2*n2*n3*n3 + n2*n2*n4*n4 + n3*n3*n4*n4) );
//Combined function (g)
//scalarvalueType gV = ( gdwV + gintV ); // Not actually needed

//Derivatives
scalarvalueType dgn1V = ( constV(2.0)*n1*(constV(1.0)-n1)*(constV(1.0)-constV(2.0)*n1) + constV(2.0)*alpha*n1*(n2*n2+n3*n3+n4*n4) );
scalarvalueType dgn2V = ( constV(2.0)*n2*(constV(1.0)-n2)*(constV(1.0)-constV(2.0)*n2) + constV(2.0)*alpha*n2*(n1*n1+n3*n3+n4*n4) );
scalarvalueType dgn3V = ( constV(2.0)*n3*(constV(1.0)-n3)*(constV(1.0)-constV(2.0)*n3) + constV(2.0)*alpha*n3*(n1*n1+n2*n2+n4*n4) );
scalarvalueType dgn4V = ( constV(2.0)*n4*(constV(1.0)-n4)*(constV(1.0)-constV(2.0)*n4) + constV(2.0)*alpha*n4*(n1*n1+n2*n2+n3*n3) );

// Residual equations
scalarvalueType rn1V = ( n1 - constV(userInputs.dtValue)*MnV*((fbV-faV)*hn1V + wV*dgn1V) );
scalarvalueType rn2V = ( n2 - constV(userInputs.dtValue)*MnV*((fbV-faV)*hn2V + wV*dgn2V) );
scalarvalueType rn3V = ( n3 - constV(userInputs.dtValue)*MnV*((fbV-faV)*hn3V + wV*dgn3V) );
scalarvalueType rn4V = ( n4 - constV(userInputs.dtValue)*MnV*((fbV-faV)*hn4V + wV*dgn4V) );
scalargradType rn1xV = ( constV(-userInputs.dtValue)*KnV*MnV*n1x );
scalargradType rn2xV = ( constV(-userInputs.dtValue)*KnV*MnV*n2x );
scalargradType rn3xV = ( constV(-userInputs.dtValue)*KnV*MnV*n3x );
scalargradType rn4xV = ( constV(-userInputs.dtValue)*KnV*MnV*n4x );
scalarvalueType rmuV = ( (constV(1.0)-hV)*facV+hV*fbcV );
scalargradType rmuxV = ( KcV*cx );
scalarvalueType rcV = ( c );
scalargradType rcxV = ( constV(-userInputs.dtValue)*McV*mux );

// Residuals for the equation to evolve the concentration 
variable_list.set_scalar_value_residual_term(0,rcV);
variable_list.set_scalar_gradient_residual_term(0,rcxV);
variable_list.set_scalar_value_residual_term(1,rmuV);
variable_list.set_scalar_gradient_residual_term(1,rmuxV);

// Residuals for the equation to evolve the order parameters 
variable_list.set_scalar_value_residual_term(2,rn1V);
variable_list.set_scalar_gradient_residual_term(2,rn1xV);
variable_list.set_scalar_value_residual_term(3,rn2V);
variable_list.set_scalar_gradient_residual_term(3,rn2xV);
variable_list.set_scalar_value_residual_term(4,rn3V);
variable_list.set_scalar_gradient_residual_term(4,rn3xV);
variable_list.set_scalar_value_residual_term(5,rn4V);
variable_list.set_scalar_gradient_residual_term(5,rn4xV);

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
