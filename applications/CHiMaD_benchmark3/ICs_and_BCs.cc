// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setInitialCondition(const dealii::Point<dim> &p,
                                            const unsigned int        index,
                                            double                   &scalar_IC,
                                            dealii::Vector<double>   &vector_IC)
{
  // ---------------------------------------------------------------------
  // ENTER THE INITIAL CONDITIONS HERE
  // ---------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the initial condition for each variable
  // according to its variable index

  // The initial condition is two circles/spheres defined
  // by a hyperbolic tangent function. The center of each circle/sphere is
  // given by "center" and its radius is given by "rad".

  double center[1][3] = {
    {0.0 / 2.0, 0.0 / 2.0, 0.0 / 2.0}
  };
  double rad[1] = {8.0};
  double dist;
  scalar_IC = 0;

  // Initial condition for the concentration field
  if (index == 0)
    {
      scalar_IC = -delta;
    }
  // Initial condition for the order parameter field
  else if (index == 1)
    {
      // Initial condition for the order parameter field
      for (unsigned int i = 0; i < 1; i++)
        {
          dist = 0.0;
          for (unsigned int dir = 0; dir < dim; dir++)
            {
              dist += (p[dir] - center[i][dir] * userInputs.domain_size[dir]) *
                      (p[dir] - center[i][dir] * userInputs.domain_size[dir]);
            }
          dist = std::sqrt(dist);

          scalar_IC += (-std::tanh((dist - rad[i]) / (0.5)));
        }
    }

  // --------------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p,
                                                  const unsigned int        index,
                                                  const unsigned int        direction,
                                                  const double              time,
                                                  double                   &scalar_BC,
                                                  dealii::Vector<double>   &vector_BC)
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
