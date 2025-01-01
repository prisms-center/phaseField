// =================================================================================
// Set the attributes of the postprocessing variables
// =================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.cc', but
// for the postprocessing expressions. It sets the attributes for each
// postprocessing expression, including its name, whether it is a vector or
// scalar (only scalars are supported at present), its dependencies on other
// variables and their derivatives, and whether to calculate an integral of the
// postprocessed quantity over the entire domain.

void
variableAttributeLoader::loadPostProcessorVariableAttributes()
{
  // Variable 0
  set_variable_name(0, "f_tot");
  set_variable_type(0, SCALAR);

  set_dependencies_value_term_RHS(0, "psi");
  set_dependencies_gradient_term_RHS(0, "");

  set_output_integral(0, true);
}

// =============================================================================================
// postProcessedFields: Set the postprocessing expressions
// =============================================================================================
// This function is analogous to 'explicitEquationRHS' and
// 'nonExplicitEquationRHS' in equations.cc. It takes in "variable_list" and
// "q_point_loc" as inputs and outputs two terms in the expression for the
// postprocessing variable -- one proportional to the test function and one
// proportional to the gradient of the test function. The index for each
// variable in this list corresponds to the index given at the top of this file
// (for submitting the terms) and the index in 'equations.cc' for assigning the
// values/derivatives of the primary variables.

template <int dim, int degree>
void
customPDE<dim, degree>::postProcessedFields(
  [[maybe_unused]] const variableContainer<dim, degree, double> &variable_list,
  [[maybe_unused]] variableContainer<dim, degree, double>       &pp_variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>     q_point_loc,
  [[maybe_unused]] const VectorizedArray<double>                 element_volume) const
{
  // --- Getting the values and derivatives of the model variables ---

  // The order parameter and its derivatives
  scalarValue psi = variable_list.get_scalar_value(4);

  // --- Setting the expressions for the terms in the postprocessing expressions
  // ---

  // --- Submitting the terms for the postprocessing expressions ---

  pp_variable_list.set_scalar_value_term_RHS(0, psi * 10e13);
}
