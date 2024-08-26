// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for
// each function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void
variableAttributeLoader::loadVariableAttributes()
{
  // Variable 0
  set_variable_name(0, "n");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "n");
  set_dependencies_gradient_term_RHS(0, "grad(n)");

  // Variable 1
  set_variable_name(1, "psi");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, TIME_INDEPENDENT);

  set_dependencies_value_term_RHS(1, "n, psi");
  set_dependencies_gradient_term_RHS(1, "grad(psi)");
  set_dependencies_value_term_LHS(1, "n, psi, change(psi)");
  set_dependencies_gradient_term_LHS(1, "grad(change(psi))");
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time
// dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a
// list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one
// proportional to the test function and one proportional to the gradient of the
// test function. The index for each variable in this list corresponds to the
// index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::explicitEquationRHS(
  variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
  dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const
{
  // --- Getting the values and derivatives of the model variables ---

  // The order parameter and its derivatives
  scalarvalueType n  = variable_list.get_scalar_value(0);
  scalargradType  nx = variable_list.get_scalar_gradient(0);

  // --- Setting the expressions for the terms in the governing equations ---

  scalarvalueType fnV   = (4.0 * n * (n - 1.0) * (n - 0.5));
  scalarvalueType eq_n  = (n - constV(userInputs.dtValue * MnV) * fnV);
  scalargradType  eqx_n = (-constV(userInputs.dtValue * KnV * MnV) * nx);

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(0, eq_n);
  variable_list.set_scalar_gradient_term_RHS(0, eqx_n);
}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time
// independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are
// not explicit time-dependent equations. It takes "variable_list" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is
// given by "q_point_loc". The function outputs two terms to variable_list --
// one proportional to the test function and one proportional to the gradient of
// the test function. The index for each variable in this list corresponds to
// the index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::nonExplicitEquationRHS(
  variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
  dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const
{
  // --- Getting the values and derivatives of the model variables ---

  scalarvalueType n = variable_list.get_scalar_value(0);

  scalarvalueType psi  = variable_list.get_scalar_value(1);
  scalargradType  psix = variable_list.get_scalar_gradient(1);

  // --- Setting the expressions for the terms in the governing equations ---

  scalarvalueType W       = constV(1.0);
  scalarvalueType p       = constV(1.5);
  scalarvalueType epsilon = constV(2.0);

  scalarvalueType eq_psi  = (W * (-psi * psi * psi + psi - 2.0 * p * psi * n * n));
  scalargradType  eqx_psi = (-epsilon * epsilon * psix);

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(1, eq_psi);
  variable_list.set_scalar_gradient_term_RHS(1, eqx_psi);
}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and
// derivatives of each of the variables at a specific quadrature point. The
// (x,y,z) location of that quadrature point is given by "q_point_loc". The
// function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function -- for the
// left-hand-side of the equation. The index for each variable in this list
// corresponds to the index given at the top of this file. If there are multiple
// elliptic equations, conditional statements should be sed to ensure that the
// correct residual is being submitted. The index of the field being solved can
// be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void
customPDE<dim, degree>::equationLHS(
  variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
  dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const
{
  // --- Getting the values and derivatives of the model variables ---

  scalarvalueType n = variable_list.get_scalar_value(0);

  scalarvalueType psi = variable_list.get_scalar_value(1);

  scalarvalueType Dpsi  = variable_list.get_change_in_scalar_value(1);
  scalargradType  Dpsix = variable_list.get_change_in_scalar_gradient(1);

  // --- Setting the expressions for the terms in the governing equations ---

  scalarvalueType W       = constV(1.0);
  scalarvalueType p       = constV(1.5);
  scalarvalueType epsilon = constV(2.0);
  scalarvalueType eq_Dpsi =
    (W * (3.0 * psi * psi * Dpsi - Dpsi + 2.0 * p * Dpsi * n * n));
  scalargradType eqx_Dpsi = (epsilon * epsilon * Dpsix);

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_LHS(1, eq_Dpsi);
  variable_list.set_scalar_gradient_term_LHS(1, eqx_Dpsi);
}
