// List of variables and residual equations for the grain growth example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0

    for (unsigned int var_index=0; var_index<10; var_index++){
        std::string var_name = "n";
        var_name.append(std::to_string(var_index+1));

        set_variable_name				(var_index,var_name);
    	set_variable_type				(var_index,SCALAR);
    	set_variable_equation_type		(var_index,PARABOLIC);

        set_dependencies_value_residual_term_RHS(var_index, "n1, n2, n3, n4, n5, n6, n7, n8, n9, n10");
        set_dependencies_gradient_residual_term_RHS(var_index, "grad(n1), grad(n2), grad(n3), grad(n4), grad(n5), grad(n6), grad(n7), grad(n8), grad(n9), grad(n10)");

    }

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
// residualNonexplicitRHS (needed only if at least one equation is elliptic)
// =================================================================================
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
