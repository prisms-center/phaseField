// List of variables and residual equations for the mechanics example
// application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void
variableAttributeLoader::loadVariableAttributes()
{
  // Variable 2
  set_variable_name(0, "u");
  set_variable_type(0, VECTOR);
  set_variable_equation_type(0, TIME_INDEPENDENT);

  set_dependencies_value_term_RHS(0, "");
  set_dependencies_gradient_term_RHS(0, "grad(u)");
  set_dependencies_value_term_LHS(0, "");
  set_dependencies_gradient_term_LHS(0, "grad(change(u))");
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
  variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
  dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc) const
{
  // --- Getting the values and derivatives of the model variables ---

  // u
  vectorgradType ux = variable_list.get_vector_gradient(0);

  // --- Setting the expressions for the terms in the governing equations ---

  vectorgradType eqx_u;

  // compute strain tensor
  dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          E[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]);
        }
    }

  // compute stress tensor
  computeStress<dim>(CIJ, E, S);

  // compute the term in the equation
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          eqx_u[i][j] = -S[i][j];
        }
    }

  // --- Submitting the terms for the governing equations ---

  variable_list.set_vector_gradient_term_RHS(0, eqx_u);
}

/// =============================================================================================
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

  // u
  vectorgradType Dux = variable_list.get_change_in_vector_gradient(0);

  // --- Setting the expressions for the terms in the governing equations ---

  vectorgradType eqx_Du;

  // compute strain tensor
  dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          E[i][j] = constV(0.5) * (Dux[i][j] + Dux[j][i]);
        }
    }

  // compute stress tensor
  computeStress<dim>(CIJ, E, S);

  // compute the term in the governing equation
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          eqx_Du[i][j] = S[i][j];
        }
    }

  // --- Submitting the terms for the governing equations ---

  variable_list.set_vector_gradient_term_LHS(0, eqx_Du);
}
