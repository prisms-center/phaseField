// List of variables and residual equations for the grain growth example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"n1");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,PARABOLIC);

	set_need_value					(0,true);
	set_need_gradient				(0,true);
	set_need_hessian				(0,false);

	set_need_value_residual_term	(0,true);
	set_need_gradient_residual_term	(0,true);

	// Variable 1
	set_variable_name				(1,"n2");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,PARABOLIC);

	set_need_value					(1,true);
	set_need_gradient				(1,true);
	set_need_hessian				(1,false);

	set_need_value_residual_term	(1,true);
	set_need_gradient_residual_term	(1,true);

	// Variable 2
	set_variable_name				(2,"n3");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,PARABOLIC);

	set_need_value					(2,true);
	set_need_gradient				(2,true);
	set_need_hessian				(2,false);

	set_need_value_residual_term	(2,true);
	set_need_gradient_residual_term	(2,true);

	// Variable 3
	set_variable_name				(3,"n4");
	set_variable_type				(3,SCALAR);
	set_variable_equation_type		(3,PARABOLIC);

	set_need_value					(3,true);
	set_need_gradient				(3,true);
	set_need_hessian				(3,false);

	set_need_value_residual_term	(3,true);
	set_need_gradient_residual_term	(3,true);

	// Variable 4
	set_variable_name				(4,"n5");
	set_variable_type				(4,SCALAR);
	set_variable_equation_type		(4,PARABOLIC);

	set_need_value					(4,true);
	set_need_gradient				(4,true);
	set_need_hessian				(4,false);

	set_need_value_residual_term	(4,true);
	set_need_gradient_residual_term	(4,true);

	// Variable 5
	set_variable_name				(5,"n6");
	set_variable_type				(5,SCALAR);
	set_variable_equation_type		(5,PARABOLIC);

	set_need_value					(5,true);
	set_need_gradient				(5,true);
	set_need_hessian				(5,false);

	set_need_value_residual_term	(5,true);
	set_need_gradient_residual_term	(5,true);

	// Variable 6
	set_variable_name				(6,"n7");
	set_variable_type				(6,SCALAR);
	set_variable_equation_type		(6,PARABOLIC);

	set_need_value					(6,true);
	set_need_gradient				(6,true);
	set_need_hessian				(6,false);

	set_need_value_residual_term	(6,true);
	set_need_gradient_residual_term	(6,true);

	// Variable 6
	set_variable_name				(7,"n8");
	set_variable_type				(7,SCALAR);
	set_variable_equation_type		(7,PARABOLIC);

	set_need_value					(7,true);
	set_need_gradient				(7,true);
	set_need_hessian				(7,false);

	set_need_value_residual_term	(7,true);
	set_need_gradient_residual_term	(7,true);

	// Variable 8
	set_variable_name				(8,"n9");
	set_variable_type				(8,SCALAR);
	set_variable_equation_type		(8,PARABOLIC);

	set_need_value					(8,true);
	set_need_gradient				(8,true);
	set_need_hessian				(8,false);

	set_need_value_residual_term	(8,true);
	set_need_gradient_residual_term	(8,true);

	// Variable 9
	set_variable_name				(9,"n10");
	set_variable_type				(9,SCALAR);
	set_variable_equation_type		(9,PARABOLIC);

	set_need_value					(9,true);
	set_need_gradient				(9,true);
	set_need_hessian				(9,false);

	set_need_value_residual_term	(9,true);
	set_need_gradient_residual_term	(9,true);

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


dealii::VectorizedArray<double> fnV = constV(0.0);
scalarvalueType ni, nj;
scalargradType nix;

// In this application, create temporary variables for the residual terms. We cannot
// call 'set_scalar_value_residual_term' and 'set_scalar_gradient_residual_term'in the
// for loop below because those functions write over the scalar value and scalar gradient
// internal variables in 'variable_list' (for performance reasons). Therefore, we wait
// to set the residual terms until all the residuals have been calculated.

std::vector<scalarvalueType> value_residuals;
value_residuals.resize(userInputs.number_of_variables);
std::vector<scalargradType> gradient_residuals;
gradient_residuals.resize(userInputs.number_of_variables);

for (unsigned int i=0; i<userInputs.number_of_variables; i++){

	ni = variable_list.get_scalar_value(i);
	nix = variable_list.get_scalar_gradient(i);
	fnV = - ni + ni*ni*ni;
	for (unsigned int j=0; j<userInputs.number_of_variables; j++){
		if (i != j){
			nj = variable_list.get_scalar_value(j);
			fnV += constV(2.0*alpha) * ni * nj*nj;
		}
	}
	value_residuals[i] = ni-constV(userInputs.dtValue*MnV)*fnV;
	gradient_residuals[i] = constV(-userInputs.dtValue*KnV*MnV)*nix;
}

for (unsigned int i=0; i<userInputs.number_of_variables; i++){
	variable_list.set_scalar_value_residual_term(i,value_residuals[i]);
	variable_list.set_scalar_gradient_residual_term(i,gradient_residuals[i]);
}

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
