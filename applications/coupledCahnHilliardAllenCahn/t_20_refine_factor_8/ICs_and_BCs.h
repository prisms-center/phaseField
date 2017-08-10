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
      
      //Constants for initial conditions equations
      double c_0 = 0.5;
      double epsilon_c = 0.05;
      double epsilon_n = 0.1;
      double psi = 1.5;

	  double dx=spanX/((double) subdivisionsX)/std::pow(2.0,refineFactor);
      double x=p[0];
      double y=p[1];

	  // Initial condition for the concentration field
	  if (index == 0){
		  if (dim == 2){
			  scalar_IC = std::cos(0.105*x)*std::cos(0.11*y);
              scalar_IC += std::cos(0.13*x)*std::cos(0.087*y)*std::cos(0.13*x)*std::cos(0.087*y);
              scalar_IC += std::cos(0.025*x-0.15*y)*std::cos(0.07*x-0.02*y);
              scalar_IC = c_0 + epsilon_c*scalar_IC;
		  }
		  else if (dim == 3) {
		  }

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
	// Third input: BC type (options are "NATURAL", "DIRICHLET", and "PERIODIC")
	// Fourth input: BC value (ignored unless the BC type is "DIRICHLET")
	// Odd inputs after the third: BC type
	// Even inputs after the third: BC value
	// Face numbering: starts at zero with the minimum of the first direction, one for the maximum of the first direction
	//						two for the minimum of the second direction, etc.

	inputBCs(0,0,"NATURAL",0);
    inputBCs(1,0,"NATURAL",0);
    inputBCs(2,0,"NATURAL",0);
    inputBCs(3,0,"NATURAL",0);
    inputBCs(4,0,"NATURAL",0);
    inputBCs(5,0,"NATURAL",0);
}


