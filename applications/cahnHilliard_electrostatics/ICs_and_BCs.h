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

  if (index == 0){
      double x=p[0];
      double y=p[1];
      double c0=0.5;
      double c1=0.04;
      
      double t1=std::cos(0.2*x)*std::cos(0.11*y);
      double t2=std::cos(0.13*x)*std::cos(0.087*y)*std::cos(0.13*x)*std::cos(0.087*y);
      double t3=std::cos(0.025*x-0.15*y)*std::cos(0.07*x-0.02*y);
      
      scalar_IC = c0 + c1*(t1+t2+t3);
  }
  else {
      scalar_IC = 0.0;
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
