// =================================================================================
// Set the attributes of the postprocessing variables
// =================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but
// for the postprocessing expressions. It sets the attributes for each
// postprocessing expression, including its name, whether it is a vector or
// scalar (only scalars are supported at present), its dependencies on other
// variables and their derivatives, and whether to calculate an integral of the
// postprocessed quantity over the entire domain.

void
variableAttributeLoader::loadPostProcessorVariableAttributes()
{
  // Variable 0
  set_variable_name(0, "mg_n");
  set_variable_type(0, SCALAR);

  set_dependencies_value_term_RHS(0, "grad(n)");
  set_dependencies_gradient_term_RHS(0, "");

  set_output_integral(0, true);

  // Variable 1
  set_variable_name(1, "f_tot");
  set_variable_type(1, SCALAR);

  set_dependencies_value_term_RHS(1, "n, grad(n)");
  set_dependencies_gradient_term_RHS(1, "");

  set_output_integral(1, true);
}

// =============================================================================================
// postProcessedFields: Set the postprocessing expressions
// =============================================================================================
// This function is analogous to 'explicitEquationRHS' and
// 'nonExplicitEquationRHS' in equations.h. It takes in "variable_list" and
// "q_point_loc" as inputs and outputs two terms in the expression for the
// postprocessing variable -- one proportional to the test function and one
// proportional to the gradient of the test function. The index for each
// variable in this list corresponds to the index given at the top of this file
// (for submitting the terms) and the index in 'equations.h' for assigning the
// values/derivatives of the primary variables.

template <int dim, int degree>
void
customPDE<dim, degree>::postProcessedFields(
  [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
    &variable_list,
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                            &pp_variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
  [[maybe_unused]] const VectorizedArray<double>             element_volume) const
{
  // --- Getting the values and derivatives of the model variables ---

  // The order parameter and its derivatives
  scalarvalueType n  = variable_list.get_scalar_value(0);
  scalargradType  nx = variable_list.get_scalar_gradient(0);

  // --- Setting the expressions for the terms in the postprocessing expressions
  // ---

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

  // --- Submitting the terms for the postprocessing expressions ---

  pp_variable_list.set_scalar_value_term_RHS(0, std::sqrt(nx[0] * nx[0] + nx[1] * nx[1]));

  pp_variable_list.set_scalar_value_term_RHS(1, f_tot);
}
