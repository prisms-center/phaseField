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

  // Initial condition parameters
  double x_denom = (1.0) * (1.0);
  double y_denom = (8.0) * (8.0);
  double z_denom = (8.0) * (8.0);

  double initial_interface_coeff = 0.08;
  double initial_radius          = 1.0;
  double c_matrix                = 1.0e-6;
  double c_precip                = 0.14;

  // set result equal to the structural order parameter initial condition
  double              r = 0.0;
  std::vector<double> ellipsoid_denoms;
  ellipsoid_denoms.push_back(x_denom);
  ellipsoid_denoms.push_back(y_denom);
  ellipsoid_denoms.push_back(z_denom);

  for (unsigned int i = 0; i < dim; i++)
    {
      r += (p(i)) * (p(i)) / ellipsoid_denoms[i];
    }
  r = sqrt(r);

  if (index == 0)
    {
      scalar_IC = 0.5 * (c_precip - c_matrix) *
                    (1.0 - std::tanh((r - initial_radius) / (initial_interface_coeff))) +
                  c_matrix;
    }
  else if (index == 1)
    {
      scalar_IC = 0.0;
    }
  else if (index == 2)
    {
      scalar_IC =
        0.5 * (1.0 - std::tanh((r - initial_radius) / (initial_interface_coeff)));
    }
  else
    {
      for (unsigned int d = 0; d < dim; d++)
        {
          vector_IC(d) = 0.0;
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
