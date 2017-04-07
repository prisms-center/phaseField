//initial condition

// =================================================================================
// Determine what (if any) files to load for the initial conditions
// =================================================================================
// Set whether PFields are enables (requires PRISMS IntegrationTools)
#define enablePFields true

// Set whether the initial conditions are loaded from file for each variable
//#define loadICs {true,true,true,false,false}
//#define loadICs {true,true,false,false,false}
#define loadICs {false,false,false,false,false}

// Set whether the file is serial or a series of parallel files
#define loadSerialFile {true,true,true,true,true}

// Set the name of the input files
#define loadFileName {"run_424a_c_1","run_424a_n1_1","run_424a_n2_1","export_test","export_test"}
//#define loadFileName {"result","result","result","export_test","export_test"}

// Set the field names in the input files
#define loadFieldName {"c","n1","n2","c","c"}

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
	  #define c_matrix 0.004 //1.0e-6
	  #define c_precip 0.14

	  double x_denom = (2.0/2.0*scaleFactor)*(2.0/2.0*scaleFactor);
	  double y_denom = (17.0/2.0*scaleFactor)*(17.0/2.0*scaleFactor);
	  double z_denom = (15.0/2.0*scaleFactor)*(15.0/2.0*scaleFactor);
	  double initial_interface_coeff = (0.02*scaleFactor);
	  double initial_radius = 1.0;

	  //set result equal to the structural order parameter initial condition
	double dx=spanX/( (double)subdivisionsX )/std::pow(2.0,refineFactor);
	double dy=spanY/( (double)subdivisionsY )/std::pow(2.0,refineFactor);
	double dz=spanZ/( (double)subdivisionsZ )/std::pow(2.0,refineFactor);
	double r=0.0;
	std::vector<double> ellipsoid_denoms;
	ellipsoid_denoms.push_back(x_denom);
	ellipsoid_denoms.push_back(y_denom);
	ellipsoid_denoms.push_back(z_denom);

	double xy_rotation = 0.0;  // 0 for n1, 2*pi/3 for n2, -2*pi/3 for n3

	double precip_center[3] = {80.0*scaleFactor, 80.0*scaleFactor, 0.0*scaleFactor};

	if (index==0){
		scalar_IC = c_matrix;
	}

	double temp1 = (p.operator()(0)-precip_center[0])*std::cos(xy_rotation) + (p.operator()(1)-precip_center[1])*std::sin(xy_rotation);
	double temp2 = (p.operator()(0)-precip_center[0])*std::sin(xy_rotation) - (p.operator()(1)-precip_center[1])*std::cos(xy_rotation);
	r = temp1*temp1/x_denom + temp2*temp2/y_denom;
	if (dim == 3){
		r += (p.operator()(2)-precip_center[2])*(p.operator()(2)-precip_center[2])/z_denom;
	}

	r = sqrt(r);

	if (index==0){
		scalar_IC += 0.5*(c_precip-c_matrix)*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));
	}

	if (index==1){
		scalar_IC += 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));
	}

	xy_rotation = 2.0*pi/3.0;  // 0 for n1, 2*pi/3 for n2, -2*pi/3 for n3

	precip_center[0] = 70.0*scaleFactor; //91.0*scaleFactor;
	precip_center[1] = 67.76*scaleFactor; //84.7*scaleFactor;
	precip_center[2] = 0.0*scaleFactor;

	temp1 = (p.operator()(0)-precip_center[0])*std::cos(xy_rotation) + (p.operator()(1)-precip_center[1])*std::sin(xy_rotation);
	temp2 = (p.operator()(0)-precip_center[0])*std::sin(xy_rotation) - (p.operator()(1)-precip_center[1])*std::cos(xy_rotation);
	r = temp1*temp1/x_denom + temp2*temp2/y_denom;
	if (dim == 3){
		r += (p.operator()(2)-precip_center[2])*(p.operator()(2)-precip_center[2])/z_denom;
	}

	r = sqrt(r);

	if (index==0){
		scalar_IC += 0.5*(c_precip-c_matrix)*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));
	}

	if (index==2){
		scalar_IC += 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));
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

	  if (index==2){
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

	this->inputBCs(0,0,"ZERO_DERIVATIVE",0);

	this->inputBCs(1,0,"ZERO_DERIVATIVE",0);

	this->inputBCs(2,0,"ZERO_DERIVATIVE",0);

	this->inputBCs(3,0,"ZERO_DERIVATIVE",0);

	if (dim == 2){
		this->inputBCs(4,0,"DIRICHLET",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0);
		this->inputBCs(4,1,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0);
	}
	else if (dim == 3){
		this->inputBCs(4,0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0);
		this->inputBCs(4,1,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0);
		this->inputBCs(4,2,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0);
	}

}


