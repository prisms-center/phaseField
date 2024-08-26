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
  set_variable_name(0, "c");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "c");
  set_dependencies_gradient_term_RHS(0, "grad(mu)");

  // Variable 1
  set_variable_name(1, "mu");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, AUXILIARY);

  set_dependencies_value_term_RHS(1, "c, phi");
  set_dependencies_gradient_term_RHS(1, "grad(c)");

  // Variable 2
  set_variable_name(2, "phi");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, TIME_INDEPENDENT);

  set_dependencies_value_term_RHS(2, "c");
  set_dependencies_gradient_term_RHS(2, "grad(phi)");
  set_dependencies_value_term_LHS(2, "");
  set_dependencies_gradient_term_LHS(2, "grad(change(phi))");
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

  // The concentration and its derivatives
  scalarvalueType c = variable_list.get_scalar_value(0);

  // The chemical potential and its derivatives
  scalargradType mux = variable_list.get_scalar_gradient(1);

  // --- Setting the expressions for the terms in the governing equations ---

  // The terms in the equations
  scalarvalueType eq_c  = c;
  scalargradType  eqx_c = constV(-McV * userInputs.dtValue) * mux;

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(0, eq_c);
  variable_list.set_scalar_gradient_term_RHS(0, eqx_c);
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

  scalarvalueType c  = variable_list.get_scalar_value(0);
  scalargradType  cx = variable_list.get_scalar_gradient(0);

  // The electric potential and its derivatives
  scalarvalueType phi  = variable_list.get_scalar_value(2);
  scalargradType  phix = variable_list.get_scalar_gradient(2);

  // --- Setting the expressions for the terms in the governing equations ---

  // The derivative of the local chemical free energy
  scalarvalueType fcV =
    2.0 * rho * (c - c_alpha) * (c_beta - c) * (c_alpha + c_beta - 2.0 * c);

  // The derivative of the local electrostatic free energy
  scalarvalueType fphiV = k * phi;

  scalarvalueType eq_mu   = fcV + fphiV;
  scalargradType  eqx_mu  = constV(KcV) * cx;
  scalarvalueType eq_phi  = -constV(-k / epsilon) * c;
  scalargradType  eqx_phi = -phix;

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(1, eq_mu);
  variable_list.set_scalar_gradient_term_RHS(1, eqx_mu);

  variable_list.set_scalar_value_term_RHS(2, eq_phi);
  variable_list.set_scalar_gradient_term_RHS(2, eqx_phi);
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
  // grad(delta phi)
  scalargradType Dphix = variable_list.get_change_in_scalar_gradient(2);

  // The residuals
  scalargradType eqx_Dphi = Dphix;

  variable_list.set_scalar_gradient_term_LHS(2, eqx_Dphi);
}
