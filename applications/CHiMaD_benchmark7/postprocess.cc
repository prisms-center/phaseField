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
  set_variable_name(0, "error_squared");
  set_variable_type(0, SCALAR);

  set_dependencies_value_term_RHS(0, "n, grad(n)");
  set_dependencies_gradient_term_RHS(0, "");

  set_output_integral(0, true);

  // Variable 1
  set_variable_name(1, "f_tot");
  set_variable_type(1, SCALAR);

  set_dependencies_value_term_RHS(1, "n, grad(n)");
  set_dependencies_gradient_term_RHS(1, "");

  set_output_integral(1, true);

  // Variable 1
  set_variable_name(2, "src");
  set_variable_type(2, SCALAR);

  set_dependencies_value_term_RHS(2, "n, grad(n)");
  set_dependencies_gradient_term_RHS(2, "");

  set_output_integral(2, true);

  // Variable 1
  set_variable_name(3, "n_sol");
  set_variable_type(3, SCALAR);

  set_dependencies_value_term_RHS(3, "n, grad(n)");
  set_dependencies_gradient_term_RHS(3, "");

  set_output_integral(3, true);
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
  const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
  variableContainer<dim, degree, dealii::VectorizedArray<double>>       &pp_variable_list,
  const dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const
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
          f_grad += constV(0.5 * kappa) * nx[i] * nx[j];
        }
    }

  // The total free energy
  f_tot = f_chem + f_grad;

  scalarvalueType source_term;
  scalarvalueType n_sol;

  scalarvalueType alpha = 0.25 + A1 * this->currentTime * std::sin(B1 * q_point_loc(0)) +
                          A2 * std::sin(B2 * q_point_loc(0) + C2 * this->currentTime);
  scalarvalueType alpha_t =
    A1 * std::sin(B1 * q_point_loc(0)) +
    A2 * C2 * std::cos(B2 * q_point_loc(0) + C2 * this->currentTime);
  scalarvalueType alpha_y =
    A1 * B1 * this->currentTime * std::cos(B1 * q_point_loc(0)) +
    A2 * B2 * std::cos(B2 * q_point_loc(0) + C2 * this->currentTime);
  scalarvalueType alpha_yy =
    -A1 * B1 * B1 * this->currentTime * std::sin(B1 * q_point_loc(0)) -
    A2 * B2 * B2 * std::sin(B2 * q_point_loc(0) + C2 * this->currentTime);

  for (unsigned i = 0; i < n.size(); i++)
    {
      source_term[i] =
        (-2.0 * std::sqrt(kappa) *
           std::tanh((q_point_loc(1)[i] - alpha[i]) / std::sqrt(2.0 * kappa)) *
           (alpha_y[i] * alpha_y[i]) +
         std::sqrt(2.0) * (alpha_t[i] - kappa * alpha_yy[i])) /
        (4.0 * std::sqrt(kappa)) /
        dealii::Utilities::fixed_power<2>(
          std::cosh((q_point_loc(1)[i] - alpha[i]) / std::sqrt(2.0 * kappa)));

      n_sol[i] =
        0.5 * (1.0 - std::tanh((q_point_loc(1)[i] - alpha[i]) / std::sqrt(2.0 * kappa)));
    }

  scalarvalueType error = (n_sol - n) * (n_sol - n);

  // --- Submitting the terms for the postprocessing expressions ---

  pp_variable_list.set_scalar_value_term_RHS(0, error);
  pp_variable_list.set_scalar_value_term_RHS(1, f_tot);
  pp_variable_list.set_scalar_value_term_RHS(2, source_term);
  pp_variable_list.set_scalar_value_term_RHS(3, n_sol);
}
