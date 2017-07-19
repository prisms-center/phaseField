template <int dim>
double InitialCondition<dim>::value (const Point<dim> &p, const unsigned int component) const
{
	  double scalar_IC = 0;
	  // =====================================================================
	  // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
	  // =====================================================================
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  // The initial condition is two circles/spheres defined
	  // by a hyperbolic tangent function. The center of each circle/sphere is
	  // given by "center" and its radius is given by "rad".

	  double center[1][3] = {{1.0/2.0,1.0/2.0,1.0/2.0}};
	  double rad[1] = {5.0};
	  double dist;
	  scalar_IC = 0;

	  // Initial condition for the concentration field
	  if (index == 0){
		  double deltaV = userInputs.get_model_constant_double("deltaV");
		  scalar_IC = -deltaV;
	  }
	  // Initial condition for the order parameter field
	  else if (index == 1) {
		  // Initial condition for the order parameter field
		  for (unsigned int i=0; i<1; i++){
			  dist = 0.0;
			  for (unsigned int dir = 0; dir < dim; dir++){
				  dist += (p[dir]-center[i][dir]*userInputs.domain_size[dir])*(p[dir]-center[i][dir]*userInputs.domain_size[dir]);
			  }
			  dist = std::sqrt(dist);

			  scalar_IC += (-std::tanh((dist-rad[i])/(0.5)));
		  }
	  }

	  // =====================================================================
	  return scalar_IC;
}

template <int dim>
void InitialConditionVec<dim>::vector_value (const dealii::Point<dim> &p, dealii::Vector<double> &vector_IC) const
{
	  // =====================================================================
	  // ENTER THE INITIAL CONDITIONS HERE FOR VECTOR FIELDS
	  // =====================================================================
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.


	  // =====================================================================
}
