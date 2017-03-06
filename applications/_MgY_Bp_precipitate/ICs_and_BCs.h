//initial condition
template <int dim>
class InitialCondition : public Function<dim>
{
public:
  unsigned int index;
  double shift;
  Vector<double> values;
  friend generalizedProblem<dim>;
  InitialCondition (const unsigned int _index, double _shift) : Function<dim>(1), index(_index), shift(_shift) {
    std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const Point<dim> &p, const unsigned int component = 0) const
  {
	  double scalar_IC = 0;
	  // =====================================================================
	  // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
	  // =====================================================================
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  // Initial condition parameters
	  #define x_denom (8.0/2.0*scaleFactor)*(8.0/2.0*scaleFactor)
	  #define y_denom (8.0/2.0*scaleFactor)*(8.0/2.0*scaleFactor)
	  #define z_denom (8.0/2.0*scaleFactor)*(8.0/2.0*scaleFactor)
	  #define initial_interface_coeff (0.02*scaleFactor)
	  #define initial_radius 1.0
	  #define c_matrix 1.0e-6
      #define c_precip 0.14
	  #define adjust_avg_c true
	  #define c_avg 0.004

	  //set result equal to the structural order parameter initial condition
	  double dx=spanX/( (double)subdivisionsX )/std::pow(2.0,refineFactor);
	  double dy=spanY/( (double)subdivisionsY )/std::pow(2.0,refineFactor);
	  double dz=spanZ/( (double)subdivisionsZ )/std::pow(2.0,refineFactor);
	  double r=0.0;
	  std::vector<double> ellipsoid_denoms;
	  ellipsoid_denoms.push_back(x_denom);
	  ellipsoid_denoms.push_back(y_denom);
	  ellipsoid_denoms.push_back(z_denom);

	  double precip_center[3] = {0.0*scaleFactor, 0.0*scaleFactor, 0.0*scaleFactor};

	  for (unsigned int i=0; i<dim; i++){
		  r += (p.operator()(i)-precip_center[i])*(p.operator()(i)-precip_center[i])/ellipsoid_denoms[i];
	  }
	  r = sqrt(r);


	  if (index==0){
		  scalar_IC = 0.5*(c_precip-c_matrix)*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff))) + c_matrix + shift;
	  }
	  else if (index==1){
		  scalar_IC = 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));
	  }

	  // =====================================================================
	  return scalar_IC;
  }
};

//initial condition
template <int dim>
class InitialConditionVec : public Function<dim>
{
public:
  unsigned int index;
  //Vector<double> values;
  InitialConditionVec (const unsigned int _index) : Function<dim>(dim), index(_index) {
    std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  void vector_value (const Point<dim> &p,Vector<double> &vector_IC) const
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
	// Third input: BC type (options are "ZERO_DERIVATIVE" and "DIRICHLET")
	// Fourth input: BC value (ignored unless the BC type is "DIRICHLET")
	// Odd inputs after the third: BC type
	// Even inputs after the third: BC value
	// Face numbering: starts at zero with the minimum of the first direction, one for the maximum of the first direction
	//						two for the minimum of the second direction, etc.

	inputBCs(0,0,"ZERO_DERIVATIVE",0);

	inputBCs(1,0,"ZERO_DERIVATIVE",0);

	if (dim == 2){
		inputBCs(2,0,"DIRICHLET",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0);
		inputBCs(2,1,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0);
	}
	else if (dim == 3){
		inputBCs(2,0,"DIRICHLET",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0);
		inputBCs(2,1,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0);
		inputBCs(2,2,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0,"ZERO_DERIVATIVE",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0,"DIRICHLET",0.0);
	}

}


