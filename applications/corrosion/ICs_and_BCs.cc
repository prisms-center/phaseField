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
  // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
  // ---------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the initial condition for each variable
  // according to its variable index

  // The initial condition is a set of overlapping circles/spheres defined
  // by a hyperbolic tangent function. The center of each circle/sphere is
  // given by "center" and its radius is given by "radius".

  double epssqV = userInputs.get_model_constant_double("epssqV");
  double deltaV = std::sqrt(2.0 * epssqV);
  double posx   = p[0];
  double posy   = p[1];
  double cx     = 0.5 * userInputs.domain_size[0];
  double cy     = userInputs.domain_size[1];
  double rad    = std::sqrt((posx - cx) * (posx - cx) + (posy - cy) * (posy - cy));
  double n0pro  = 0.5 * (1.0 - std::tanh((rad0 - rad) / deltaV));

  if (index == 0)
    {
      scalar_IC = n0pro;
      if (scalar_IC > 1.0)
        scalar_IC = 1.0;
    }
  else if (index == 2)
    {
      scalar_IC = 1.0 - n0pro;
      if (scalar_IC > 1.0)
        scalar_IC = 1.0;
    }
  else if (index == 5)
    {
      scalar_IC = 1000.0;
    }
  else if (index == 6)
    {
      scalar_IC = 0.0;
    }
  else
    {
      scalar_IC = 0.0;
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
  /*
if (index == 2){
  if (direction == 1){
      double x=p[0];
      double y=p[1];
      scalar_BC=std::sin(y/7.0);
  }
}
  */
  // -------------------------------------------------------------------------
}
