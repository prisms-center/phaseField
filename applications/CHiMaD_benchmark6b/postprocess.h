// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){

	// Variable 0
	set_variable_name				(0,"f_tot");
	set_variable_type				(0,SCALAR);

	set_need_value_residual_term	(0,true);
	set_need_gradient_residual_term	(0,false);

    set_output_integral         	(0,true);

}

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

	scalarvalueType f_tot = constV(0.0);

	// The concentration and its derivatives 
	scalarvalueType c = variable_list.get_scalar_value(0);
	scalargradType cx = variable_list.get_scalar_gradient(0);
	scalarvalueType phi = variable_list.get_scalar_value(2);
    
	// The homogenous chemical free energy
	scalarvalueType f_chem = rho*(c-c_alpha)*(c-c_alpha)*(c_beta-c)*(c_beta-c);

    // The homogenous electrostaric free energy
    scalarvalueType f_elec = 0.5*k*c*phi;

	// The gradient free energy
	scalarvalueType f_grad = constV(0.0);

	for (int i=0; i<dim; i++){
		for (int j=0; j<dim; j++){
			f_grad += constV(0.5*KcV)*cx[i]*cx[j];
		}
	}

	// The total free energy
	f_tot = f_chem + f_elec + f_grad;


// Residuals for the equation to evolve the order parameter 
pp_variable_list.set_scalar_value_residual_term(0, f_tot);


}
