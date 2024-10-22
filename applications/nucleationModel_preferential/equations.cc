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
  set_dependencies_gradient_term_RHS(0, "c, grad(c), n, grad(n)");

  set_allowed_to_nucleate(0, false);
  set_need_value_nucleation(0, true);

  // Variable 1
  set_variable_name(1, "n");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(1, "c, n");
  set_dependencies_gradient_term_RHS(1, "grad(n)");

  set_allowed_to_nucleate(1, true);
  set_need_value_nucleation(1, true);
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
  scalarvalueType c  = variable_list.get_scalar_value(0);
  scalargradType  cx = variable_list.get_scalar_gradient(0);

  // The order parameter and its derivatives
  scalarvalueType n  = variable_list.get_scalar_value(1);
  scalargradType  nx = variable_list.get_scalar_gradient(1);

  // --- Setting the expressions for the terms in the governing equations ---

  // Interpolation function and its derivative
  scalarvalueType hV  = 3.0 * n * n - 2.0 * n * n * n;
  scalarvalueType hnV = 6.0 * n - 6.0 * n * n;

  // KKS model c_alpha and c_beta as a function of c and h
  scalarvalueType c_alpha =
    (B2 * (c - cbtmin * hV) + A2 * calmin * hV) / (A2 * hV + B2 * (1.0 - hV));
  scalarvalueType c_beta = (A2 * (c - calmin * (1.0 - hV)) + B2 * cbtmin * (1.0 - hV)) /
                           (A2 * hV + B2 * (1.0 - hV));

  // Free energy for each phase and their first and second derivatives
  scalarvalueType faV  = A0 + A2 * (c_alpha - calmin) * (c_alpha - calmin);
  scalarvalueType fbV  = B0 + B2 * (c_beta - cbtmin) * (c_beta - cbtmin);
  scalarvalueType fbcV = 2.0 * B2 * (c_beta - cbtmin);

  // Double-Well function (can be used to tune the interfacial energy)
  scalarvalueType fbarriernV = 2.0 * n - 6.0 * n * n + 4.0 * n * n * n;

  // -------------------------------------------------
  // Nucleation expressions
  // -------------------------------------------------
  VectorizedArray<double> source_term = constV(0.0);
  VectorizedArray<double> gamma       = constV(1.0);
  seedNucleus(q_point_loc, source_term, gamma);
  // -------------------------------------------------

  // Set the terms in the governing equations

  // For concentration
  scalarvalueType eq_c = c;
  scalargradType  eqx_c =
    constV(-McV * userInputs.dtValue) * (cx + (c_alpha - c_beta) * hnV * nx);
  // For order parameter (gamma is a variable order parameter mobility factor)
  scalarvalueType eq_n =
    n - constV(userInputs.dtValue * MnV) * gamma *
          ((fbV - faV) * hnV - (c_beta - c_alpha) * fbcV * hnV + W_barrier * fbarriernV);
  scalargradType eqx_n = constV(-userInputs.dtValue * KnV * MnV) * gamma * nx;

  // --- Submitting the terms for the governing equations ---

  // Terms for the equation to evolve the concentration
  variable_list.set_scalar_value_term_RHS(0, eq_c);
  variable_list.set_scalar_gradient_term_RHS(0, eqx_c);

  // Terms for the equation to evolve the order parameter
  variable_list.set_scalar_value_term_RHS(1, eq_n + source_term);
  variable_list.set_scalar_gradient_term_RHS(1, eqx_n);
}

// =================================================================================
// seedNucleus: a function particular to this app
// =================================================================================
template <int dim, int degree>
void
customPDE<dim, degree>::seedNucleus(
  const Point<dim, VectorizedArray<double>> &q_point_loc,
  VectorizedArray<double>                   &source_term,
  VectorizedArray<double>                   &gamma) const
{
  // Loop through all of the seeded nuclei
  for (const auto &thisNucleus : this->nuclei)
    {
      if (thisNucleus.seededTime + thisNucleus.seedingTime > this->currentTime)
        {
          // Calculate the weighted distance function to the order parameter
          // freeze boundary (weighted_dist = 1.0 on that boundary)
          VectorizedArray<double> weighted_dist = this->weightedDistanceFromNucleusCenter(
            thisNucleus.center,
            userInputs.get_nucleus_freeze_semiaxes(thisNucleus.orderParameterIndex),
            q_point_loc,
            thisNucleus.orderParameterIndex);

          for (unsigned i = 0; i < gamma.size(); i++)
            {
              if (weighted_dist[i] <= 1.0)
                {
                  gamma[i] = 0.0;

                  // Seed a nucleus if it was added to the list of nuclei this
                  // time step
                  if (thisNucleus.seedingTimestep == this->currentIncrement)
                    {
                      // Find the weighted distance to the outer edge of the
                      // nucleus and use it to calculate the order parameter
                      // source term (r = 1.0 on that boundary)
                      Point<dim, double> q_point_loc_element;
                      for (unsigned int j = 0; j < dim; j++)
                        {
                          q_point_loc_element(j) = q_point_loc(j)[i];
                        }
                      double r = this->weightedDistanceFromNucleusCenter(
                        thisNucleus.center,
                        userInputs.get_nucleus_semiaxes(thisNucleus.orderParameterIndex),
                        q_point_loc_element,
                        thisNucleus.orderParameterIndex);

                      double avg_semiaxis = 0.0;
                      for (unsigned int j = 0; j < dim; j++)
                        {
                          avg_semiaxis += thisNucleus.semiaxes[j];
                        }
                      avg_semiaxis /= dim;

                      source_term[i] =
                        0.5 *
                        (1.0 - std::tanh(avg_semiaxis * (r - 1.0) / interface_coeff));
                    }
                }
            }
        }
    }
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
{}

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
