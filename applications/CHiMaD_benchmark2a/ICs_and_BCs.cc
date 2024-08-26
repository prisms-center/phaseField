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

  // Constants for initial conditions equations
  double c_0       = 0.5;
  double epsilon_c = 0.05;
  double epsilon_n = 0.1;
  double psi       = 1.5;

  double dx = userInputs.domain_size[0] / ((double) userInputs.subdivisions[0]) /
              std::pow(2.0, userInputs.refine_factor);
  double x = p[0];
  double y = p[1];

  // Initial condition for the concentration field
  if (index == 0)
    {
      if (dim == 2)
        {
          scalar_IC = std::cos(0.105 * x) * std::cos(0.11 * y);
          scalar_IC += std::cos(0.13 * x) * std::cos(0.087 * y) * std::cos(0.13 * x) *
                       std::cos(0.087 * y);
          scalar_IC += std::cos(0.025 * x - 0.15 * y) * std::cos(0.07 * x - 0.02 * y);
          scalar_IC = c_0 + epsilon_c * scalar_IC;
        }
      else if (dim == 3)
        {
        }
    }
  // Initial condition for the chemical potential field
  if (index == 1)
    {
      if (dim == 2)
        {
          scalar_IC = 0.0;
        }
      else if (dim == 3)
        {
        }
    }
  // Initial condition for order parameters
  if (index >= 2)
    {
      double j = ((double) index) - 1.0;
      if (dim == 2)
        {
          double term1;
          double term2;
          term1 = std::cos(0.01 * j * x - 4.0) * std::cos((0.007 + 0.01 * j) * y);
          term1 += std::cos((0.11 + 0.01 * j) * x) * std::cos((0.11 + 0.01 * j) * y);
          term2 = std::cos((0.046 + 0.001 * j) * x + (0.0405 + 0.001 * j) * y) *
                  std::cos((0.031 + 0.001 * j) * x - (0.004 + 0.001 * j) * y);
          term2     = psi * term2 * term2;
          scalar_IC = term1 + term2;
          scalar_IC = epsilon_n * (scalar_IC * scalar_IC);
        }
      else if (dim == 3)
        {
        }
    }

  // ---------------------------------------------------------------------
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
