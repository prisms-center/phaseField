//initial condition

// =================================================================================

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

	  // Initial condition parameters
	  double pi = (2.0*std::acos(0.0));

	  #define adjust_avg_c false
	  #define c_avg 0.004
	  #define c_matrix 0.02 //1.0e-6
	  #define c_precip 0.125
	  // Initial condition parameters
	#define x_denom (2.5)*(2.5)
	#define y_denom (2.5)*(2.5)
	#define z_denom (2.5)*(2.5)
	#define initial_interface_coeff interface_coeff
	#define initial_radius 1.0

	  //set result equal to the structural order parameter initial condition
	  double r=0.0;
	  std::vector<double> ellipsoid_denoms;
	  ellipsoid_denoms.push_back(x_denom);
	  ellipsoid_denoms.push_back(y_denom);
	  ellipsoid_denoms.push_back(z_denom);

	  double precip_center[3] = {spanX/2.0, spanY/2.0, spanZ/2.0};

	  for (unsigned int i=0; i<dim; i++){
		  r += (p.operator()(i)-precip_center[i])*(p.operator()(i)-precip_center[i])/ellipsoid_denoms[i];
	  }
	  r = sqrt(r);



	if (index==0){
		scalar_IC = c_matrix;
		//scalar_IC = 0.5*(c_precip-c_matrix)*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff))) + c_matrix;
	}
	else if (index == 1){
		scalar_IC = 0.0;
		//scalar_IC = 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));

	}
	else if (index==2||index==3){
		scalar_IC = 0.0;
	}


	  // =====================================================================
	  return scalar_IC;
  }


//initial condition
template <int dim>
  void InitialConditionVec<dim>::vector_value (const Point<dim> &p,Vector<double> &vector_IC) const
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

template <int dim,int degree>
void customPDE<dim,degree>::setBCs(){
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

	this->inputBCs(0,0,"PERIODIC",0);

	this->inputBCs(1,0,"PERIODIC",0);

	this->inputBCs(2,0,"PERIODIC",0);

	this->inputBCs(3,0,"PERIODIC",0);

	if (dim == 2){
		this->inputBCs(4,0,"PERIODIC",0.0);
		this->inputBCs(4,1,"PERIODIC",0.0);
	}
	else if (dim == 3){
		this->inputBCs(4,0,"PERIODIC",0.0);
		this->inputBCs(4,1,"PERIODIC",0.0);
		this->inputBCs(4,2,"PERIODIC",0.0);
	}

}


