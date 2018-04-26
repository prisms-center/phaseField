// ===========================================================================
// FUNCTIONS FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim>
double InitialCondition<dim>::value (const dealii::Point<dim> &p, const unsigned int component) const
{
  double scalar_IC;
  // ---------------------------------------------------------------------
  // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
  // ---------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the initial condition for each variable
  // according to its variable index.


  if (index ==0){
      return 0.5*(1.0-std::tanh((p[0]-50.0)/1.5));
  }
  else {
      return p[0]/200.0;
  }

}

template <int dim>
void InitialConditionVec<dim>::vector_value (const dealii::Point<dim> &p, dealii::Vector<double> &vector_IC) const
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
double NonUniformDirichletBC<dim>::value (const dealii::Point<dim> &p, const unsigned int component) const
{
    double scalar_BC=0;
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
void NonUniformDirichletBCVec<dim>::vector_value (const dealii::Point<dim> &p, dealii::Vector<double> &vector_BC) const
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
