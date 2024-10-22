// =============================================================================================
// loadPostProcessorVariableAttributes: Set the attributes of the postprocessing
// variables
// =============================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but
// for the postprocessing expressions. It sets the attributes for each
// postprocessing expression, including its name, whether it is a vector or
// scalar (only scalars are supported at present), its dependencies on other
// variables and their derivatives, and whether to calculate an integral of the
// postprocessed quantity over the entire domain. Note: this function is not a
// member of customPDE.

void
variableAttributeLoader::loadPostProcessorVariableAttributes()
{
  // Variable 0
  set_variable_name(0, "f_tot");
  set_variable_type(0, SCALAR);

  set_dependencies_value_term_RHS(0, "c, grad(c), n, grad(n), biharm");
  set_dependencies_gradient_term_RHS(0, "");

  set_output_integral(0, true);
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
  scalarvalueType c = variable_list.get_scalar_value(0);

  // n1
  scalarvalueType n  = variable_list.get_scalar_value(1);
  scalargradType  nx = variable_list.get_scalar_gradient(1);

  // biharm
  scalarvalueType biharm = variable_list.get_scalar_value(2);

  // --- Setting the expressions for the terms in the postprocessing expressions
  // ---

  scalarvalueType f_tot = constV(0.0);

  scalarvalueType faV = 0.5 * c * c / 16.0;
  scalarvalueType fbV = 0.5 * (c - 1.0) * (c - 1.0) / 16.0;
  scalarvalueType hV  = 3.0 * n * n - 2.0 * n * n * n;

  scalarvalueType normgradn = std::sqrt(nx.norm_square());
  scalargradType  normal    = nx / (normgradn + constV(1.0e-16));

  scalarvalueType gamma;
  scalargradType  dgammadnormal;
  anisotropy(normal, gamma, dgammadnormal);

  scalarvalueType f_chem = (constV(1.0) - hV) * faV + hV * fbV;

  // anisotropy code
  scalarvalueType f_grad = constV(0.5) * gamma * gamma * nx * nx;

  scalarvalueType f_reg = constV(0.5 * delta2) * biharm * biharm;

  f_tot = f_chem + f_grad + f_reg;

  // end anisotropy code
  f_tot = f_chem + f_grad;

  // --- Submitting the terms for the postprocessing expressions ---
  pp_variable_list.set_scalar_value_term_RHS(0, f_tot);
}
