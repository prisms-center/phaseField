
template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// The order parameter and its derivatives (names here should match those in the macros above)
scalargradType nx = variable_list.get_scalar_gradient(0);

scalargradType pp_field;
pp_field[0] = nx[0];
pp_field[1] = nx[1];


// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
pp_variable_list.set_scalar_value_residual_term(0, std::sqrt(pp_field[0]*pp_field[0]+pp_field[1]*pp_field[1]));


}
