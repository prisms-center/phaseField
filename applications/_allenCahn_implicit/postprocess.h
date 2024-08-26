// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void
variableAttributeLoader::loadPostProcessorVariableAttributes()
{
  set_variable_name(0, "f_tot");
  set_variable_type(0, SCALAR);

  set_output_integral(0, true);

  set_dependencies_value_residual_term_RHS(0, "psi, grad(psi)");
  set_dependencies_gradient_residual_term_RHS(0, "");
}

template <int dim, int degree>
void
customPDE<dim, degree>::postProcessedFields(
  const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
  variableContainer<dim, degree, dealii::VectorizedArray<double>>       &pp_variable_list,
  const dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const
{
  // The order parameter and its derivatives
  scalarvalueType n  = variable_list.get_scalar_value(0);
  scalargradType  nx = variable_list.get_scalar_gradient(0);

  scalargradType pp_field;
  pp_field[0] = nx[0];
  pp_field[1] = nx[1];

  scalarvalueType f_tot = constV(0.0);

  // The homogenous free energy
  scalarvalueType f_chem = (n * n * n * n - 2.0 * n * n * n + n * n);

  // The gradient free energy
  scalarvalueType f_grad = constV(0.0);

  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += constV(0.5 * KnV) * nx[i] * nx[j];
        }
    }

  // The total free energy
  f_tot = f_chem + f_grad;

  // Residuals for the equation to evolve the order parameter
  pp_variable_list.set_scalar_value_residual_term(0, f_tot);
}
