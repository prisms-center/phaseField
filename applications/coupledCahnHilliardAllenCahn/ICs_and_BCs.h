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
      
      //Constants for initial conditions equations
      double c_0 = 0.5;
      double epsilon_c = 0.05;
      double epsilon_n = 0.1;
      double psi = 1.5;

<<<<<<< HEAD
	  double dx=spanX/((double) subdivisionsX)/std::pow(2.0,refineFactor);
      double x=p[0];
      double y=p[1];
=======
	  // The initial condition is two circles/spheres defined
	  // by a hyperbolic tangent function. The center of each circle/sphere is
	  // given by "center" and its radius is given by "rad".

	  double center[2][3] = {{1.0/3.0,1.0/3.0,1.0/3.0},{3.0/4.0,3.0/4.0,3.0/4.0}};
	  double rad[12] = {userInputs.domain_size[0]/5.0, userInputs.domain_size[0]/12.0};
	  double dist;
	  scalar_IC = 0;
>>>>>>> 0147129c14b77ac7320a8032fc57a37558811c95

	  if (index == 0){
<<<<<<< HEAD
		  if (dim == 2){
			  scalar_IC = std::cos(0.105*x)*std::cos(0.11*y);
              scalar_IC += std::cos(0.13*x)*std::cos(0.087*y)*std::cos(0.13*x)*std::cos(0.087*y);
              scalar_IC += std::cos(0.025*x-0.15*y)*std::cos(0.07*x-0.02*y);
              scalar_IC = c_0 + epsilon_c*scalar_IC;
		  }
		  else if (dim == 3) {
=======
		  scalar_IC = 0.009;
	  }

	  for (unsigned int i=0; i<2; i++){
		  dist = 0.0;
		  for (unsigned int dir = 0; dir < dim; dir++){
			  dist += (p[dir]-center[i][dir]*userInputs.domain_size[dir])*(p[dir]-center[i][dir]*userInputs.domain_size[dir]);
>>>>>>> 0147129c14b77ac7320a8032fc57a37558811c95
		  }
		  dist = std::sqrt(dist);

<<<<<<< HEAD
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
      
=======
		  // Initial condition for the concentration field
		  if (index == 0){
			  scalar_IC += 0.5*(0.125)*(1.0-std::tanh((dist-rad[i])/(1.0)));
		  }
		  else {
			  scalar_IC += 0.5*(1.0-std::tanh((dist-rad[i])/(1.0)));
		  }
	  }

>>>>>>> 0147129c14b77ac7320a8032fc57a37558811c95
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
<<<<<<< HEAD
  }
};

template <int dim>
void generalizedProblem<dim>::setBCs(){

	// =====================================================================
	// ENTER THE BOUNDARY CONDITIONS HERE
	// =====================================================================
	// This function sets the BCs for the problem variables
	// The function "inputBCs" should be called for each component of
	// each variable and should be in numerical order. Four input arguments
	// set the same BC on the entire boundary. Two plus two times the
	// number of dimensions inputs sets separate BCs on each face of the domain.
	// Inputs to "inputBCs":
	// First input: variable number
	// Second input: component number
	// Third input: BC type (options are "ZERO_DERIVATIVE", "DIRICHLET", and "PERIODIC")
	// Fourth input: BC value (ignored unless the BC type is "DIRICHLET")
	// Odd inputs after the third: BC type
	// Even inputs after the third: BC value
	// Face numbering: starts at zero with the minimum of the first direction, one for the maximum of the first direction
	//						two for the minimum of the second direction, etc.

	inputBCs(0,0,"ZERO_DERIVATIVE",0);
    inputBCs(1,0,"ZERO_DERIVATIVE",0);
    inputBCs(2,0,"ZERO_DERIVATIVE",0);
    inputBCs(3,0,"ZERO_DERIVATIVE",0);
    inputBCs(4,0,"ZERO_DERIVATIVE",0);
    inputBCs(5,0,"ZERO_DERIVATIVE",0);
=======
>>>>>>> 0147129c14b77ac7320a8032fc57a37558811c95
}
