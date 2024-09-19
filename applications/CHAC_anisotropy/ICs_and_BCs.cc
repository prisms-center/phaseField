// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                                            [[maybe_unused]] const unsigned int index,
                                            [[maybe_unused]] double            &scalar_IC,
                                            [[maybe_unused]] Vector<double>    &vector_IC)
{
  // ---------------------------------------------------------------------
  // ENTER THE INITIAL CONDITIONS HERE
  // ---------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the initial condition for each variable
  // according to its variable index

  double r = 0.0;

  if (index == 0)
    {
      r = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++)
        {
          r += (p[dir] - userInputs.domain_size[dir] / 2.0) *
               (p[dir] - userInputs.domain_size[dir] / 2.0);
        }
      r         = std::sqrt(r);
      double n  = 0.5 * (1.0 - std::tanh((r - userInputs.domain_size[0] / 4.0) / 4.0));
      scalar_IC = 0.082 * 16.0 / (userInputs.domain_size[0] / 4.0) +
                  (3.0 * n * n - 2.0 * n * n * n);
    }
  else
    {
      r = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++)
        {
          r += (p[dir] - userInputs.domain_size[dir] / 2.0) *
               (p[dir] - userInputs.domain_size[dir] / 2.0);
        }
      r         = std::sqrt(r);
      scalar_IC = 0.5 * (1.0 - std::tanh((r - userInputs.domain_size[0] / 4.0) / 4.0));
    }

  // --------------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setNonUniformDirichletBCs(
  [[maybe_unused]] const Point<dim>  &p,
  [[maybe_unused]] const unsigned int index,
  [[maybe_unused]] const unsigned int direction,
  [[maybe_unused]] const double       time,
  [[maybe_unused]] double            &scalar_BC,
  [[maybe_unused]] Vector<double>    &vector_BC)
{
  // --------------------------------------------------------------------------
  // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
  // --------------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the boundary condition for each variable
  // according to its variable index. This function can be left blank if there
  // are no non-uniform Dirichlet boundary conditions. For BCs that change in
  // time, you can access the current time through the variable "time". The
  // boundary index can be accessed via the variable "direction", which starts
  // at zero and uses the same order as the BC specification in parameters.in
  // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).

  // -------------------------------------------------------------------------
}
