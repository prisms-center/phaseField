// ===========================================================================
// FUNCTIONS FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim>
double
InitialCondition<dim>::value(const dealii::Point<dim> &p,
                             const unsigned int        component) const
{
  double scalar_IC;
  // ---------------------------------------------------------------------
  // ENTER THE INITIAL CONDITIONS HERE
  // ---------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the initial condition for each variable
  // according to its variable index.

  // The initial condition is a set of overlapping circles/spheres defined
  // by a hyperbolic tangent function. The center of each circle/sphere is
  // given by "center" and its radius is given by "radius".

  double center[12][3] = {
    {0.1, 0.3,  0},
    {0.8, 0.7,  0},
    {0.5, 0.2,  0},
    {0.4, 0.4,  0},
    {0.3, 0.9,  0},
    {0.8, 0.1,  0},
    {0.9, 0.5,  0},
    {0.0, 0.1,  0},
    {0.1, 0.6,  0},
    {0.5, 0.6,  0},
    {1,   1,    0},
    {0.7, 0.95, 0}
  };
  double rad[12] = {12, 14, 19, 16, 11, 12, 17, 15, 20, 10, 11, 14};
  double dist;
  scalar_IC = 0;
  for (unsigned int i = 0; i < 12; i++)
    {
      dist = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++)
        {
          dist += (p[dir] - center[i][dir] * userInputs.domain_size[dir]) *
                  (p[dir] - center[i][dir] * userInputs.domain_size[dir]);
        }
      dist = std::sqrt(dist);

      scalar_IC += 0.5 * (1.0 - std::tanh((dist - rad[i]) / 1.5));
    }
  if (scalar_IC > 1.0)
    scalar_IC = 1.0;

  // ---------------------------------------------------------------------
  return scalar_IC;
}

template <int dim>
void
InitialConditionVec<dim>::vector_value(const dealii::Point<dim> &p,
                                       dealii::Vector<double>   &vector_IC) const
{
  // ---------------------------------------------------------------------
  // ENTER THE INITIAL CONDITIONS HERE FOR VECTOR FIELDS
  // ---------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the initial condition for each variable
  // according to its variable index.

  // ---------------------------------------------------------------------
}

// ===========================================================================
// FUNCTIONS FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim>
double
NonUniformDirichletBC<dim>::value(const dealii::Point<dim> &p,
                                  const unsigned int        component) const
{
  double scalar_BC = 0;
  // --------------------------------------------------------------------------
  // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE FOR SCALAR FIELDS
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
  return scalar_BC;
}

template <int dim>
void
NonUniformDirichletBCVec<dim>::vector_value(const dealii::Point<dim> &p,
                                            dealii::Vector<double>   &vector_BC) const
{
  // --------------------------------------------------------------------------
  // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE FOR VECTOR FIELDS
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
