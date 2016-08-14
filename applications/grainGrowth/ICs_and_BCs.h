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
	  double scalar_IC;
	  // =====================================================================
	  // ENTER THE INITIAL CONDITIONS HERE FOR SCALAR FIELDS
	  // =====================================================================
	  // Enter the function describing conditions for the fields at point "p".
	  // Use "if" statements to set the initial condition for each variable
	  // according to its variable index.

	  //return  0.5+ 0.2*(0.5 - (double)(std::rand() % 100 )/100.0);

	  double dx=spanX/((double) subdivisionsX)/std::pow(2.0,refineFactor);
	  double x_loc[10] = {0.3, 0.9, 0.5, 0.4, 0.3, 0.0, 0.0, 0.6, 0.7, 1.0};
	  double y_loc[10] = {0.3, 0.7, 0.1, 0.5, 0.9, 0.0, 0.7, 0.6, 0.4, 0.1};
	  double rad[10] =   {6, 7, 9, 8, 5, 6, 4, 3, 8, 5};
	  double dist;
	  scalar_IC = 0;

	  #if problemDIM == 2
	  	  dist = p.distance(Point<dim>(x_loc[index]*spanX,y_loc[index]*spanY));
	  #elif problemDIM == 3
	  	  dist = p.distance(Point<dim>(x_loc[index]*spanX,y_loc[index]*spanY,0.5*spanZ));
	  #endif

	  if (dist < rad[index]){
		  scalar_IC = 1.0;
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
	// Third input: BC type (options are "ZERO_DERIVATIVE" and "DIRICHLET")
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
	inputBCs(6,0,"ZERO_DERIVATIVE",0);
	inputBCs(7,0,"ZERO_DERIVATIVE",0);
	inputBCs(8,0,"ZERO_DERIVATIVE",0);
	inputBCs(9,0,"ZERO_DERIVATIVE",0);

	// =====================================================================

}


