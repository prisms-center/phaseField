template <int dim>
class InitialCondition : public Function<dim>
{
public:
  unsigned int index;
  Vector<double> values;
  InitialCondition (const unsigned int _index) : Function<dim>(1), index(_index) {
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
    
    if (index == 0){
      double epsilon=0.01;
      double c0 = 0.5;
      double dx=spanX/((double) subdivisionsX)/std::pow(2.0,refineFactor);
      scalar_IC = c0+epsilon*(std::cos(0.105*p[0])*std::cos(0.11*p[1])+
	std::pow(std::cos(0.13*p[0])*std::cos(0.087*p[1]),2.0)+
	std::cos(0.025*p[0]-0.15*p[1])*std::cos(0.07*p[0]-0.02*p[1]));

    }
    else {
      scalar_IC = 0.0;
    }
    
    // =====================================================================
    return scalar_IC;
  }
};

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
	// Third input: BC type (options are "ZERO_DERIVATIVE", "DIRICHLET", and "PERIODIC")
	// Fourth input: BC value (ignored unless the BC type is "DIRICHLET")
	// Odd inputs after the third: BC type
	// Even inputs after the third: BC value
	// Face numbering: starts at zero with the minimum of the first direction, one for the maximum of the first direction
	//						two for the minimum of the second direction, etc.
  
  inputBCs(0,0,"ZERO_DERIVATIVE",0);
  inputBCs(1,0,"ZERO_DERIVATIVE",0);
  
}
