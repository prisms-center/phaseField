// ===========================================================================
// FUNCTIONS FOR INITIAL CONDITIONS
// ===========================================================================


template <int dim>
double InitialCondition<dim>::value (const Point<dim> &p, const unsigned int component) const
{
	  double scalar_IC = 0;
	  /// --------------------------------------------------------------------------
	  // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
	  // --------------------------------------------------------------------------
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  double dx=userInputs.domain_size[0]/((double) userInputs.subdivisions[0])/std::pow(2.0,userInputs.refine_factor);
	  double r=0.0;

	  // Initial condition for the concentration field
	  if (index == 0){
		  if (dim == 2){
              scalar_IC = userInputs.get_model_constant_double("c_avg");;
		  }
		  else if (dim == 3) {
              scalar_IC = userInputs.get_model_constant_double("c_avg");;
		  }

	  }
	  // Initial condition for the structural order parameter field
	  else {
		  if (dim == 2){
              scalar_IC = 0.0;
		  }
		  else if (dim == 3){
              scalar_IC = 0.0;
		  }
	  }

	  // --------------------------------------------------------------------------
	  return scalar_IC;
}

template <int dim>
void InitialConditionVec<dim>::vector_value (const dealii::Point<dim> &p, dealii::Vector<double> &vector_IC) const
{
	  // --------------------------------------------------------------------------
	  // ENTER THE INITIAL CONDITIONS HERE FOR VECTOR FIELDS
	  // --------------------------------------------------------------------------
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.


	  // --------------------------------------------------------------------------
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
    // are no non-uniform Dirichlet boundary conditions.


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
    // are no non-uniform Dirichlet boundary conditions.


    // -------------------------------------------------------------------------

}
