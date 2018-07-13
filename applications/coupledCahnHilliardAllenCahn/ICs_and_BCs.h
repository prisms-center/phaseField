// ===========================================================================
// FUNCTIONS FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim>
double InitialCondition<dim>::value (const Point<dim> &p, const unsigned int component) const
{
	  double scalar_IC = 0;
	  // --------------------------------------------------------------------------
	  // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
	  // --------------------------------------------------------------------------
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  // The initial condition is two circles/spheres defined
	  // by a hyperbolic tangent function. The center of each circle/sphere is
	  // given by "center" and its radius is given by "rad".

      dealii::Tensor<1,dim> center1 = userInputs.get_model_constant_rank_1_tensor("center1");
      dealii::Tensor<1,dim> center2 = userInputs.get_model_constant_rank_1_tensor("center2");
      double radius1 = userInputs.get_model_constant_double("radius1");
      double radius2 = userInputs.get_model_constant_double("radius2");
      double matrix_concentration = userInputs.get_model_constant_double("matrix_concentration");



	  double dist;
	  scalar_IC = 0.0;

      dist = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++){
          dist += (p[dir]-center1[dir])*(p[dir]-center1[dir]);
      }
      dist = std::sqrt(dist);

      // Initial condition for the concentration field
      if (index == 0){
          scalar_IC += matrix_concentration + (1.0-matrix_concentration)*0.5*(1.0-std::tanh((dist-radius1)/(1.0)));
      }
      else {
          scalar_IC += 0.5*(1.0-std::tanh((dist-radius1)/(1.0)));
      }

      dist = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++){
          dist += (p[dir]-center2[dir])*(p[dir]-center2[dir]);
      }
      dist = std::sqrt(dist);

      // Initial condition for the concentration field
      if (index == 0){
          scalar_IC += (1.0-matrix_concentration)*0.5*(1.0-std::tanh((dist-radius2)/(1.0)));
      }
      else {
          scalar_IC += 0.5*(1.0-std::tanh((dist-radius2)/(1.0)));
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
