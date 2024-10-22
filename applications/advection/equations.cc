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
  set_variable_equation_type(0, IMPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "n, old(n), grad(n)");
  set_dependencies_value_term_LHS(0, "change(n), grad(change(n))");
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
  [[maybe_unused]] VectorizedArray<double> element_volume) const
{}

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
  [[maybe_unused]] VectorizedArray<double> element_volume) const
{
  // Getting necessary variables
  scalarvalueType n     = variable_list.get_scalar_value(0);
  scalargradType  nx    = variable_list.get_scalar_gradient(0);
  scalarvalueType n_old = variable_list.get_old_scalar_value(0);

  // Converting prescibed velocity to a vectorvalueType
  vectorvalueType vel;
  for (unsigned int i = 0; i < dim; i++)
    {
      vel[i] = constV(velocity[i]);
    }

  // Norm of the local velocity
  scalarvalueType u_l2norm = 1.0e-12 + vel.norm_square();

  // Submission terms
  scalarvalueType residual = (n_old - n - constV(userInputs.dtValue) * vel * nx);
  scalarvalueType eq_n     = residual;

  variable_list.set_scalar_value_term_RHS(0, eq_n);
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
  [[maybe_unused]] VectorizedArray<double> element_volume) const
{
  // Getting necessary variables
  scalarvalueType change_n  = variable_list.get_change_in_scalar_value(0);
  scalargradType  change_nx = variable_list.get_change_in_scalar_gradient(0);

  // Converting prescibed velocity to a vectorvalueType
  vectorvalueType vel;
  for (unsigned int i = 0; i < dim; i++)
    {
      vel[i] = constV(velocity[i]);
    }

  // Norm of the local velocity
  scalarvalueType u_l2norm = 1.0e-12 + vel.norm_square();

  // Submission terms
  scalarvalueType residual = (change_n + constV(userInputs.dtValue) * vel * change_nx);
  scalarvalueType eq_n     = residual;

  variable_list.set_scalar_value_term_LHS(0, eq_n);
}
