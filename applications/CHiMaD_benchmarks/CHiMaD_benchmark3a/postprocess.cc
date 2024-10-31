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
  // Variable 1
  set_variable_name(0, "f_tot");
  set_variable_type(0, SCALAR);

  set_dependencies_value_term_RHS(0, "u, phi, grad(phi)");
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

  // The temperature and its derivatives
  scalarvalueType u = variable_list.get_scalar_value(0);

  // The order parameter and its derivatives
  scalarvalueType phi  = variable_list.get_scalar_value(1);
  scalargradType  phix = variable_list.get_scalar_gradient(1);

  // --- Setting the expressions for the terms in the postprocessing expressions
  // ---

  double lambda = (D / 0.6267 / W0 / W0);

  scalarvalueType f_tot = constV(0.0);

  // The homogenous free energy
  scalarvalueType f_chem =
    -0.5 * phi * phi + 0.25 * phi * phi * phi * phi +
    lambda * u * phi * (1.0 - 2.0 / 3.0 * phi * phi + 1.0 / 5.0 / phi * phi * phi * phi);

  // The azimuthal angle
  scalarvalueType theta;
  for (unsigned i = 0; i < phi.size(); i++)
    {
      theta[i] = std::atan2(phix[1][i], phix[0][i]);
    }

  scalarvalueType W =
    constV(W0) *
    (constV(1.0) + constV(epsilonM) * std::cos(constV(mult) * (theta - constV(theta0))));

  // The gradient free energy
  scalarvalueType f_grad = constV(0.5) * W * W * (phix[0] * phix[0] + phix[1] * phix[1]);

  // The total free energy
  f_tot = f_chem + f_grad;

  // --- Submitting the terms for the postprocessing expressions ---
  pp_variable_list.set_scalar_value_term_RHS(0, f_tot);
}
