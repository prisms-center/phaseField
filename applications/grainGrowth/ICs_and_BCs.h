template <int dim>
double InitialCondition<dim>::value (const dealii::Point<dim> &p, const unsigned int component) const
{
  double scalar_IC;
	  // =====================================================================
	  // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
	  // =====================================================================
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  // The initial condition is a set of overlapping circles/spheres defined
	  // by a hyperbolic tangent function. The center of each circle/sphere is
	  // given by "center" and its radius is given by "radius".

	  double center[10][3] = {{0.3,0.3,0},{0.7,0.7,0},{0.5,0.1,0},{0.4,0.5,0},{0.3,0.9,0},{0.1,0.1,0},{0.1,0.7,0},{0.6,0.6,0},{0.7,0.4,0},{0.7,0.2,0}};      
	  double rad[12] = {6, 7, 9, 8, 5, 6, 4, 3, 8, 5};
	  double dist = 0.0;
	  scalar_IC = 0;

	  for (unsigned int dir = 0; dir < dim; dir++){
		  dist += (p[dir]-center[index][dir]*userInputs.domain_size[dir])*(p[dir]-center[index][dir]*userInputs.domain_size[dir]);
	  }
	  dist = std::sqrt(dist);

	  scalar_IC +=	0.5*(1.0-std::tanh((dist-rad[index])/0.8));

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
