// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){

	// Variable 0
	set_variable_name				(0,"diff1");
	set_variable_type				(0,SCALAR);

	set_need_value_residual_term	(0,true);
	set_need_gradient_residual_term	(0,true);

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
