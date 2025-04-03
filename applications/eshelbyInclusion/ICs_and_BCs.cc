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

  for (unsigned int d = 0; d < dim; d++)
    {
      vector_IC(d) = 0.0;
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
  //radius is declared in equations.cc but this code doesn't see that so I'm defining it here as well
  for (unsigned int i = 0; i < dim; i++)
  {
    //double radius = 10.0; //radius of inclusion
    double mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]); //magnitude function
    double A = (incRadius*incRadius*incRadius)/(6*0.9); //All constants for the displacement equation
    //eigenstrain is only a function of distance from inclusion
    //double sfts = 0.01 * (0.5 + 0.5 * (1.0 - std::exp(-20.0 * (mag - incRadius))) /
    //     (1.0 + std::exp(-20.0 * (mag - incRadius))));
    double sfts = 0.01 * (0.5 + 0.5 * (-1.0 * std::tanh(-20.0 * (mag - incRadius))));
    if (direction == 1)
    {
      //eigenstrain is uniform and diagonal so I need only three directional functions
      //directional functions change based on direction of displacement
      double xdir = 0.8*(p[0]/mag) + 3*(p[0]*p[0]*p[0])/(mag*mag*mag); // with a 0.1 poisson
      double ydir = 0.8*(-1*p[0]/mag) + 3*(p[0]*p[1]*p[1])/(mag*mag*mag);
      double zdir = 0.8*(-1*p[0]/mag) + 3*(p[0]*p[2]*p[2])/(mag*mag*mag);
      vector_BC(i) = A*(1/(mag*mag))*sfts*(xdir + ydir + zdir);
    }
    if (direction == 3)
    {
      double xdir = 0.8*(-1*p[1]/mag) + 3*(p[1]*p[0]*p[0])/(mag*mag*mag);
      double ydir = 0.8*(p[1]/mag) + 3*(p[1]*p[1]*p[1])/(mag*mag*mag);
      double zdir = 0.8*(-1*p[1]/mag) + 3*(p[1]*p[2]*p[2])/(mag*mag*mag);
      vector_BC(i) = A*(1/(mag*mag))*sfts*(xdir + ydir + zdir);
    }
    if (direction == 5)
    {
      double xdir = 0.8*(-1*p[2]/mag) + 3*(p[2]*p[0]*p[0])/(mag*mag*mag);
      double ydir = 0.8*(-1*p[2]/mag) + 3*(p[2]*p[1]*p[1])/(mag*mag*mag);
      double zdir = 0.8*(p[2]/mag) + 3*(p[2]*p[2]*p[2])/(mag*mag*mag);
     vector_BC(i) = A*(1/(mag*mag))*sfts*(xdir + ydir + zdir);
    }
  }
}
