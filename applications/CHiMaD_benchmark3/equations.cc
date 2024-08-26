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
  set_variable_name(0, "u");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "u,mu,grad(phi)");
  set_dependencies_gradient_term_RHS(0, "grad(u)");

  // Variable 1
  set_variable_name(1, "phi");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(1, "phi,mu");
  set_dependencies_gradient_term_RHS(1, "");

  // Variable 2
  set_variable_name(2, "mu");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, AUXILIARY);

  set_dependencies_value_term_RHS(2, "phi,u,grad(phi)");
  set_dependencies_gradient_term_RHS(2, "grad(phi)");
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

  // The temperature and its derivatives
  scalarvalueType u  = variable_list.get_scalar_value(0);
  scalargradType  ux = variable_list.get_scalar_gradient(0);

  // The order parameter and its derivatives
  scalarvalueType phi  = variable_list.get_scalar_value(1);
  scalargradType  phix = variable_list.get_scalar_gradient(1);

  // The order parameter chemical potential and its derivatives
  scalarvalueType mu = variable_list.get_scalar_value(2);

  // --- Setting the expressions for the terms in the governing equations ---

  // The azimuthal angle
  scalarvalueType theta;
  for (unsigned i = 0; i < phi.size(); i++)
    {
      theta[i] = std::atan2(phix[1][i], phix[0][i]);
    }

  // Anisotropic gradient energy coefficient, its derivative and square
  scalarvalueType W =
    constV(W0) *
    (constV(1.0) + constV(epsilonM) * std::cos(constV(mult) * (theta - constV(theta0))));
  scalarvalueType tau = W / constV(W0);

  // Define terms in the equations
  scalarvalueType eq_u   = (u + constV(0.5) * mu * constV(userInputs.dtValue) / tau);
  scalargradType  eqx_u  = (constV(-D * userInputs.dtValue) * ux);
  scalarvalueType eq_phi = (phi + constV(userInputs.dtValue) * mu / tau);

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(0, eq_u);
  variable_list.set_scalar_gradient_term_RHS(0, eqx_u);

  variable_list.set_scalar_value_term_RHS(1, eq_phi);
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

  // The temperature and its derivatives
  scalarvalueType u = variable_list.get_scalar_value(0);

  // The order parameter and its derivatives
  scalarvalueType phi  = variable_list.get_scalar_value(1);
  scalargradType  phix = variable_list.get_scalar_gradient(1);

  // --- Setting the expressions for the terms in the governing equations ---

  // The coupling constant, determined from solvability theory
  double lambda = (D / 0.6267 / W0 / W0);

  // Derivative of the free energy density with respect to phi
  scalarvalueType f_phi =
    -(phi - constV(lambda) * u * (constV(1.0) - phi * phi)) * (constV(1.0) - phi * phi);

  // The azimuthal angle
  scalarvalueType theta;
  for (unsigned i = 0; i < phi.size(); i++)
    {
      theta[i] = std::atan2(phix[1][i], phix[0][i]);
    }

  // Anisotropic gradient energy coefficient, its derivative and square
  scalarvalueType W =
    constV(W0) *
    (constV(1.0) + constV(epsilonM) * std::cos(constV(mult) * (theta - constV(theta0))));
  scalarvalueType W_theta =
    constV(-W0) *
    (constV(epsilonM) * constV(mult) * std::sin(constV(mult) * (theta - constV(theta0))));
  scalarvalueType tau = W / constV(W0);

  // The anisotropy term that enters in to the  equation for mu
  scalargradType aniso;
  aniso[0] = W * W * phix[0] - W * W_theta * phix[1];
  aniso[1] = W * W * phix[1] + W * W_theta * phix[0];

  // Define the terms in the equations
  scalarvalueType eq_mu  = (-f_phi);
  scalargradType  eqx_mu = (-aniso);

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(2, eq_mu);
  variable_list.set_scalar_gradient_term_RHS(2, eqx_mu);
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
{}
