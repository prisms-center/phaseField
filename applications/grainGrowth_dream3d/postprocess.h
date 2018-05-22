// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){

	// Variable 0
	set_variable_name				(0,"diff1");
	set_variable_type				(0,SCALAR);

    // set_dependencies_value_residual_term_RHS(0, "n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14");
    // set_dependencies_gradient_residual_term_RHS(0, "grad(n0), grad(n1), grad(n2), grad(n3), grad(n4), grad(n5), grad(n6), grad(n7), grad(n8), grad(n9), grad(n10), grad(n11), grad(n12), grad(n13), grad(n14)");

    set_dependencies_value_residual_term_RHS(0, "n0, n1, n2, n3, n4, n5");
    set_dependencies_gradient_residual_term_RHS(0, "grad(n0), grad(n1), grad(n2), grad(n3), grad(n4), grad(n5)");

    set_output_integral         	(0,true);

}

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

dealii::VectorizedArray<double> fnV = constV(0.0);
scalarvalueType ni, nj;
scalargradType nix;

for (unsigned int i=0; i<userInputs.number_of_variables; i++){
	ni = variable_list.get_scalar_value(i) ;
	nix = variable_list.get_scalar_gradient(i);
	fnV = - ni + ni*ni*ni;
	for (unsigned int j=0; j<userInputs.number_of_variables; j++){
		if (i != j){
			nj = variable_list.get_scalar_value(j);
			fnV += constV(2.0*alpha) * ni * nj*nj;
		}
	}
	if (i==1){
		pp_variable_list.set_scalar_value_residual_term(0, -constV(userInputs.dtValue*MnV)*fnV);
		pp_variable_list.set_scalar_gradient_residual_term(0, constV(-userInputs.dtValue*KnV*MnV)*nix);
	}
}
// Residuals for the equation to evolve the order parameter



}
