
// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"mg_c");
	set_variable_type				(0,SCALAR);

	set_need_value_residual_term	(0,true);
	set_need_gradient_residual_term	(0,false);

    set_output_integral         	(0,true);

}

// =================================================================================
// Define the expressions for the post-processed fields
// =================================================================================

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// The order parameter and its derivatives 
scalargradType cx = variable_list.get_scalar_gradient(0);

scalargradType pp_field;
pp_field[0] = cx[0];
pp_field[1] = cx[1];


// Residuals for the equation to evolve the order parameter 
pp_variable_list.set_scalar_value_residual_term(0, std::sqrt(pp_field[0]*pp_field[0]+pp_field[1]*pp_field[1]));

}
