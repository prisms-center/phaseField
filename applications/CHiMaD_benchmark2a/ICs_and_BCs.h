// ===========================================================================
// FUNCTIONS FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim>
double InitialCondition<dim>::value (const Point<dim> &p, const unsigned int component) const
{
	  double scalar_IC = 0;

	  // ---------------------------------------------------------------------
      // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
      // ---------------------------------------------------------------------
      // Enter the function describing conditions for the fields at point "p".
      // Use "if" statements to set the initial condition for each variable
      // according to its variable index.

      //Constants for initial conditions equations
      double c_0 = 0.5;
      double epsilon_c = 0.05;
      double epsilon_n = 0.1;
      double psi = 1.5;

	  double dx=userInputs.domain_size[0]/((double) userInputs.subdivisions[0])/std::pow(2.0,userInputs.refine_factor);
      double x=p[0];
      double y=p[1];

	  // Initial condition for the concentration field
	  if (index == 0){
		  if (dim == 2){
			  scalar_IC = std::cos(0.105*x)*std::cos(0.11*y);
			  scalar_IC += std::cos(0.13*x)*std::cos(0.087*y)*std::cos(0.13*x)*std::cos(0.087*y);
			  scalar_IC += std::cos(0.025*x-0.15*y)*std::cos(0.07*x-0.02*y);
			  scalar_IC = c_0 + epsilon_c*scalar_IC;
		  }
		  else if (dim == 3) {
		  }

	  }
	  // Initial condition for the chemical potential field
	  if (index == 1){
		  if (dim == 2){
			  scalar_IC = 0.0;
		  }
		  else if (dim == 3) {
		  }
	  }
	  // Initial condition for order parameters
	  if (index >= 2){
		  double j = ((double) index)-1.0;
		  if (dim == 2){
			  double term1;
			  double term2;
			  term1 = std::cos(0.01*j*x-4.0)*std::cos((0.007+0.01*j)*y);
			  term1 += std::cos((0.11+0.01*j)*x)*std::cos((0.11+0.01*j)*y);
			  term2 = std::cos((0.046+0.001*j)*x + (0.0405+0.001*j)*y)*std::cos((0.031+0.001*j)*x - (0.004+0.001*j)*y);
			  term2 = psi*term2*term2;
			  scalar_IC = term1 + term2;
			  scalar_IC = epsilon_n*(scalar_IC*scalar_IC);
		  }
		  else if (dim == 3) {
		  }
	  }

	  // ---------------------------------------------------------------------
	  return scalar_IC;
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
