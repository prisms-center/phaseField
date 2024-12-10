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
  set_variable_name(0, "f_tot");
  set_variable_type(0, SCALAR);

  set_dependencies_value_term_RHS(0, "c,n,grad(n)");
  set_dependencies_gradient_term_RHS(0, "");

  set_output_integral(0, true);

  // Variable 0
  set_variable_name(1, "c_grad");
  set_variable_type(1, SCALAR);

  set_dependencies_value_term_RHS(1, "grad(c)");
  set_dependencies_gradient_term_RHS(1, "");

  set_output_integral(1, false);
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

  // c
  scalarvalueType c  = variable_list.get_scalar_value(0);
  scalargradType  cx = variable_list.get_scalar_gradient(0);

  // n
  scalarvalueType n  = variable_list.get_scalar_value(1);
  scalargradType  nx = variable_list.get_scalar_gradient(1);

  // --- Setting the expressions for the terms in the postprocessing expressions
  // ---

  // Free energy for each phase and their first and second derivatives
  scalarvalueType fa =
    (-1.6704 - 4.776 * c + 5.1622 * c * c - 2.7375 * c * c * c + 1.3687 * c * c * c * c);
  scalarvalueType fb = (5.0 * c * c - 5.9746 * c - 1.5924);

  // Interpolation function and its derivative
  scalarvalueType h = (10.0 * n * n * n - 15.0 * n * n * n * n + 6.0 * n * n * n * n * n);

  // The homogenous free energy
  scalarvalueType f_chem = (constV(1.0) - h) * fa + h * fb;

  // The gradient free energy
  scalarvalueType f_grad = constV(0.5 * Kn) * nx * nx;

  // The total free energy
  scalarvalueType f_tot;
  f_tot = f_chem + f_grad;

  // The magnitude of the gradient of c
  scalarvalueType mag_grad_c = constV(0.0);
  for (unsigned int i = 0; i < dim; i++)
    {
      mag_grad_c = mag_grad_c + cx[i] * cx[i];
    }
  for (unsigned int v = 0; v < c.size(); v++)
    {
      mag_grad_c[v] = sqrt(mag_grad_c[v]);
    }

  // --- Submitting the terms for the postprocessing expressions ---

  pp_variable_list.set_scalar_value_term_RHS(0, f_tot);
  pp_variable_list.set_scalar_value_term_RHS(1, mag_grad_c);
}
