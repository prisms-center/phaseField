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
  // For the input file 'parameters.in'
  for (unsigned int var_index = 0; var_index < 8; var_index++)
    {
      // For the input file 'parameters_large_2D.in'
      // for (unsigned int var_index=0; var_index<12; var_index++){
      std::string var_name = "n";
      var_name.append(std::to_string(var_index));

      set_variable_name(var_index, var_name);
      set_variable_type(var_index, SCALAR);
      set_variable_equation_type(var_index, EXPLICIT_TIME_DEPENDENT);

      // For the input file 'parameters.in'
      set_dependencies_value_term_RHS(var_index, "n0, n1, n2, n3, n4, n5, n6, n7");
      set_dependencies_gradient_term_RHS(
        var_index,
        "grad(n0), grad(n1), grad(n2), grad(n3), grad(n4), grad(n5), grad(n6), "
        "grad(n7)");

      // For the input file 'parameters_large_2D.in'
      // set_dependencies_value_term_RHS(var_index, "n0, n1, n2, n3, n4, n5, n6,
      // n7 , n8, n9, n10, n11"); set_dependencies_gradient_term_RHS(var_index,
      // "grad(n0), grad(n1), grad(n2), grad(n3), grad(n4), grad(n5), grad(n6),
      // grad(n7), grad(n8), grad(n9), grad(n10), grad(n11)");
    }
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

  VectorizedArray<double> fnV = constV(0.0);
  scalarvalueType         ni, nj;
  scalargradType          nix;

  // In this application, create temporary variables for the residual terms. We
  // cannot call 'set_scalar_value_residual_term' and
  // 'set_scalar_gradient_residual_term'in the for loop below because those
  // functions write over the scalar value and scalar gradient internal
  // variables in 'variable_list' (for performance reasons). Therefore, we wait
  // to set the residual terms until all the residuals have been calculated.

  std::vector<scalarvalueType> value_terms;
  value_terms.resize(userInputs.number_of_variables);
  std::vector<scalargradType> gradient_terms;
  gradient_terms.resize(userInputs.number_of_variables);

  for (unsigned int i = 0; i < userInputs.number_of_variables; i++)
    {
      ni  = variable_list.get_scalar_value(i);
      nix = variable_list.get_scalar_gradient(i);
      fnV = -ni + ni * ni * ni;
      for (unsigned int j = 0; j < userInputs.number_of_variables; j++)
        {
          if (i != j)
            {
              nj = variable_list.get_scalar_value(j);
              fnV += constV(2.0 * alpha) * ni * nj * nj;
            }
        }
      value_terms[i]    = ni - constV(userInputs.dtValue * MnV) * fnV;
      gradient_terms[i] = constV(-userInputs.dtValue * KnV * MnV) * nix;
    }

  // --- Submitting the terms for the governing equations ---

  for (unsigned int i = 0; i < userInputs.number_of_variables; i++)
    {
      variable_list.set_scalar_value_term_RHS(i, value_terms[i]);
      variable_list.set_scalar_gradient_term_RHS(i, gradient_terms[i]);
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
