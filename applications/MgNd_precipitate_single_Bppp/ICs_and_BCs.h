// ===========================================================================
// FUNCTIONS FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim>
double InitialCondition<dim>::value (const dealii::Point<dim> &p, const unsigned int component) const
{
  double scalar_IC=0;
      // --------------------------------------------------------------------------
      // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
      // --------------------------------------------------------------------------
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  // Initial condition parameters
      double x_denom = (1.0)*(1.0);
      double y_denom = (8.0)*(8.0);
      double z_denom = (8.0)*(8.0);

	  double initial_interface_coeff = 0.08;
	  double initial_radius = 1.0;
	  double c_matrix = 1.0e-6;
      double c_precip = 0.14;

	  //set result equal to the structural order parameter initial condition
	  double r=0.0;
	  std::vector<double> ellipsoid_denoms;
	  ellipsoid_denoms.push_back(x_denom);
	  ellipsoid_denoms.push_back(y_denom);
	  ellipsoid_denoms.push_back(z_denom);

	  for (unsigned int i=0; i<dim; i++){
		  r += (p(i))*(p(i))/ellipsoid_denoms[i];
	  }
	  r = sqrt(r);


	  if (index==0){
		  scalar_IC = 0.5*(c_precip-c_matrix)*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff))) + c_matrix;
	  }
      else if (index==1){
          scalar_IC = 0.0;
      }
	  else if (index==2){
		  scalar_IC = 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));
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

	  if (index==3){
		  vector_IC(0) = 0.0;
		  vector_IC(1) = 0.0;
          if (dim == 3){
        	  vector_IC(2) = 0.0;
          }
	  }
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
