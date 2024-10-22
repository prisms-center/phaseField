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

  set_dependencies_value_term_RHS(1, "c, n1, n2, n3, n4");
  set_dependencies_gradient_term_RHS(1, "grad(c)");

  // Variable 2
  set_variable_name(2, "n1");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(2, "c, n1, n2, n3, n4");
  set_dependencies_gradient_term_RHS(2, "grad(n1)");

  // Variable 3
  set_variable_name(3, "n2");
  set_variable_type(3, SCALAR);
  set_variable_equation_type(3, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(3, "c, n1, n2, n3, n4");
  set_dependencies_gradient_term_RHS(3, "grad(n2)");

  // Variable 4
  set_variable_name(4, "n3");
  set_variable_type(4, SCALAR);
  set_variable_equation_type(4, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(4, "c, n1, n2, n3, n4");
  set_dependencies_gradient_term_RHS(4, "grad(n3)");

  // Variable 5
  set_variable_name(5, "n4");
  set_variable_type(5, SCALAR);
  set_variable_equation_type(5, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(5, "c, n1, n2, n3, n4");
  set_dependencies_gradient_term_RHS(5, "grad(n4)");
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
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // --- Getting the values and derivatives of the model variables ---

  // The concentration and its derivatives
  scalarvalueType c = variable_list.get_scalar_value(0);

  // The chemical potential and its derivatives
  scalargradType mux = variable_list.get_scalar_gradient(1);

  // The order parameters and their derivatives
  scalarvalueType n1  = variable_list.get_scalar_value(2);
  scalargradType  n1x = variable_list.get_scalar_gradient(2);
  scalarvalueType n2  = variable_list.get_scalar_value(3);
  scalargradType  n2x = variable_list.get_scalar_gradient(3);
  scalarvalueType n3  = variable_list.get_scalar_value(4);
  scalargradType  n3x = variable_list.get_scalar_gradient(4);
  scalarvalueType n4  = variable_list.get_scalar_value(5);
  scalargradType  n4x = variable_list.get_scalar_gradient(5);

  // --- Setting the expressions for the terms in the governing equations ---

  // Free energy for each phase and their first and second derivatives
  scalarvalueType faV = (constV(2.0) * (c - constV(0.3)) * (c - constV(0.3)));
  scalarvalueType fbV = (constV(2.0) * (c - constV(0.7)) * (c - constV(0.7)));

  // Interpolation function and its derivatives
  scalarvalueType hn1V =
    (n1 * n1 * (constV(30.0) * n1 * n1 - constV(60.0) * n1 + constV(30.0)));
  scalarvalueType hn2V =
    (n2 * n2 * (constV(30.0) * n2 * n2 - constV(60.0) * n2 + constV(30.0)));
  scalarvalueType hn3V =
    (n3 * n3 * (constV(30.0) * n3 * n3 - constV(60.0) * n3 + constV(30.0)));
  scalarvalueType hn4V =
    (n4 * n4 * (constV(30.0) * n4 * n4 - constV(60.0) * n4 + constV(30.0)));

  // Derivatives
  scalarvalueType dgn1V =
    (constV(2.0) * n1 * (constV(1.0) - n1) * (constV(1.0) - constV(2.0) * n1) +
     constV(2.0) * alpha * n1 * (n2 * n2 + n3 * n3 + n4 * n4));
  scalarvalueType dgn2V =
    (constV(2.0) * n2 * (constV(1.0) - n2) * (constV(1.0) - constV(2.0) * n2) +
     constV(2.0) * alpha * n2 * (n1 * n1 + n3 * n3 + n4 * n4));
  scalarvalueType dgn3V =
    (constV(2.0) * n3 * (constV(1.0) - n3) * (constV(1.0) - constV(2.0) * n3) +
     constV(2.0) * alpha * n3 * (n1 * n1 + n2 * n2 + n4 * n4));
  scalarvalueType dgn4V =
    (constV(2.0) * n4 * (constV(1.0) - n4) * (constV(1.0) - constV(2.0) * n4) +
     constV(2.0) * alpha * n4 * (n1 * n1 + n2 * n2 + n3 * n3));

  // Terms for the governing equations
  scalarvalueType eq_n1 =
    (n1 - constV(userInputs.dtValue) * MnV * ((fbV - faV) * hn1V + wV * dgn1V));
  scalarvalueType eq_n2 =
    (n2 - constV(userInputs.dtValue) * MnV * ((fbV - faV) * hn2V + wV * dgn2V));
  scalarvalueType eq_n3 =
    (n3 - constV(userInputs.dtValue) * MnV * ((fbV - faV) * hn3V + wV * dgn3V));
  scalarvalueType eq_n4 =
    (n4 - constV(userInputs.dtValue) * MnV * ((fbV - faV) * hn4V + wV * dgn4V));
  scalargradType eqx_n1 = (constV(-userInputs.dtValue) * KnV * MnV * n1x);
  scalargradType eqx_n2 = (constV(-userInputs.dtValue) * KnV * MnV * n2x);
  scalargradType eqx_n3 = (constV(-userInputs.dtValue) * KnV * MnV * n3x);
  scalargradType eqx_n4 = (constV(-userInputs.dtValue) * KnV * MnV * n4x);

  scalarvalueType eq_c  = (c);
  scalargradType  eqx_c = (constV(-userInputs.dtValue) * McV * mux);

  // --- Submitting the terms for the governing equations ---

  // Terms for the equations to evolve the concentration
  variable_list.set_scalar_value_term_RHS(0, eq_c);
  variable_list.set_scalar_gradient_term_RHS(0, eqx_c);

  // Terms for the equations to evolve the order parameters
  variable_list.set_scalar_value_term_RHS(2, eq_n1);
  variable_list.set_scalar_gradient_term_RHS(2, eqx_n1);
  variable_list.set_scalar_value_term_RHS(3, eq_n2);
  variable_list.set_scalar_gradient_term_RHS(3, eqx_n2);
  variable_list.set_scalar_value_term_RHS(4, eq_n3);
  variable_list.set_scalar_gradient_term_RHS(4, eqx_n3);
  variable_list.set_scalar_value_term_RHS(5, eq_n4);
  variable_list.set_scalar_gradient_term_RHS(5, eqx_n4);
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
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // --- Getting the values and derivatives of the model variables ---

  scalarvalueType c  = variable_list.get_scalar_value(0);
  scalargradType  cx = variable_list.get_scalar_gradient(0);

  // The order parameters and their derivatives
  scalarvalueType n1 = variable_list.get_scalar_value(2);
  scalarvalueType n2 = variable_list.get_scalar_value(3);
  scalarvalueType n3 = variable_list.get_scalar_value(4);
  scalarvalueType n4 = variable_list.get_scalar_value(5);

  // --- Setting the expressions for the terms in the governing equations ---

  // Free energy for each phase and their first and second derivatives
  scalarvalueType facV = (constV(4.0) * (c - constV(0.3)));
  scalarvalueType fbcV = (constV(4.0) * (c - constV(0.7)));

  // Interpolation function and its derivatives
  scalarvalueType hV =
    (n1 * n1 * n1 * (constV(6.0) * n1 * n1 - constV(15.0) * n1 + constV(10.0)) +
     n2 * n2 * n2 * (constV(6.0) * n2 * n2 - constV(15.0) * n2 + constV(10.0)) +
     n3 * n3 * n3 * (constV(6.0) * n3 * n3 - constV(15.0) * n3 + constV(10.0)) +
     n4 * n4 * n4 * (constV(6.0) * n4 * n4 - constV(15.0) * n4 + constV(10.0)));

  // The terms for the governing equations
  scalarvalueType eq_mu  = ((constV(1.0) - hV) * facV + hV * fbcV);
  scalargradType  eqx_mu = (KcV * cx);

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(1, eq_mu);
  variable_list.set_scalar_gradient_term_RHS(1, eqx_mu);
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
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{}
