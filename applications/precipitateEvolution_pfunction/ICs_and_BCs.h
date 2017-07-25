//initial condition
template <int dim>
double InitialCondition<dim>::value (const dealii::Point<dim> &p, const unsigned int component) const
{
  double scalar_IC=0.0;
	  // =====================================================================
	  // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
	  // =====================================================================
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  double center[4][3] = {{1.0/3.0,1.0/3.0,0.5},{2.0/3.0,2.0/3.0,0.5},{3.0/4.0,1.0/4.0,0.5},{1.0/4.0,3.0/4,0.5}};
	  double rad[4] = {userInputs.domain_size[0]/16.0, userInputs.domain_size[0]/16.0, userInputs.domain_size[0]/16.0, userInputs.domain_size[0]/16.0};
	  double orientation[4] = {1,1,2,3};
	  double dx=userInputs.domain_size[0]/((double) userInputs.subdivisions[0])/std::pow(2.0,userInputs.refine_factor);
	  double dist;
	  scalar_IC = 0;

	  if (index==0){
		  scalar_IC = 0.04;
	  }

	  for (unsigned int i=0; i<4; i++){
		  dist = 0.0;
		  for (unsigned int dir = 0; dir < dim; dir++){
			  dist += (p[dir]-center[i][dir]*userInputs.domain_size[dir])*(p[dir]-center[i][dir]*userInputs.domain_size[dir]);
		  }
		  dist = std::sqrt(dist);

		  if (index == orientation[i]){
			  scalar_IC +=	0.5*(1.0-std::tanh((dist-rad[i])/(dx)));
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

	  if (index==4){
		  vector_IC(0) = 0.0;
		  vector_IC(1) = 0.0;
          if (dim == 3){
        	  vector_IC(2) = 0.0;
          }
	  }
	  // =====================================================================
}
