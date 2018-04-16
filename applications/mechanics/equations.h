// List of variables and residual equations for the mechanics example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){

	// Variable 2
	set_variable_name				(0,"u");
	set_variable_type				(0,VECTOR);
	set_variable_equation_type		(0,ELLIPTIC);

	set_need_value					(0,false);
	set_need_gradient				(0,true);
	set_need_hessian				(0,false);

	set_need_value_residual_term	(0,false);
	set_need_gradient_residual_term	(0,true);

	set_need_value_LHS				(0,false);
	set_need_gradient_LHS			(0,false);
	set_need_hessian_LHS			(0,false);


    set_need_value_change_LHS		(0,false);
	set_need_gradient_change_LHS	(0,true);
	set_need_hessian_change_LHS		(0,false);

    set_need_value_residual_term_LHS	(0,false);
	set_need_gradient_residual_term_LHS	(0,true);

    set_equations_are_nonlinear(false);

}


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
template <int dim, int degree>
void customPDE<dim,degree>::residualRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//u
vectorgradType ux = variable_list.get_vector_gradient(0);
vectorgradType ruxV;

//compute strain tensor
dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		E[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]);
	}
}

//compute stress tensor
computeStress<dim>(CIJ, E, S);

//compute residual
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		ruxV[i][j] = -S[i][j];
	}
}

variable_list.set_vector_gradient_residual_term(0,ruxV);

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
template <int dim, int degree>
void customPDE<dim,degree>::residualLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {


//u
vectorgradType Dux = variable_list.get_change_in_vector_gradient(0);
vectorgradType ruxV;

//compute strain tensor
dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		E[i][j]= constV(0.5)*(Dux[i][j]+Dux[j][i]);
	}
}

//compute stress tensor
computeStress<dim>(CIJ, E, S);

//compute residual
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		ruxV[i][j] = S[i][j];
	}
}

variable_list.set_vector_gradient_residual_term_LHS(0,ruxV);

}
