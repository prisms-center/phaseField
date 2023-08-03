// Unit test(s) for the method "computeInvM"
template <int dim, int degree>
class testInvM: public MatrixFreePDE<dim,degree>
{
 public:
  testInvM(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) {

	  // Initialize the test field object (needed for computeInvM())
	  Field<dim> test_field(SCALAR,EXPLICIT_TIME_DEPENDENT,"c");
	  this->fields.push_back(test_field);

	  //init the MatrixFreePDE class for testing
	  this->initForTests();

	  //call computeInvM()
	  this->computeInvM();
	  invMNorm=this->invM.l2_norm();

  };
  ~testInvM(){
      this->matrixFreeObject.clear();
  };

  double invMNorm;

  void setBCs(){};

  // Function to set the initial conditions (in ICs_and_BCs.h)
  void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC){};

  // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
  void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC){};

 private:
  //RHS implementation for explicit solve
  void getRHS(const MatrixFree<dim,double> &data,
	      std::vector<vectorType*> &dst,
	      const std::vector<vectorType*> &src,
	      const std::pair<unsigned int,unsigned int> &cell_range) const{};

  // Function to set the RHS of the governing equations for explicit time dependent equations (in equations.h)
  void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                   dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {};

  // Function to set the RHS of the governing equations for all other equations (in equations.h)
  void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                   dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {};

  // Function to set the LHS of the governing equations (in equations.h)
  void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                   dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {};


};

template <int dim,typename T>
  bool unitTest<dim,T>::test_computeInvM(int argc, char** argv, userInputParameters<dim> userInputs){
  	bool pass = false;
	std::cout << "\nTesting 'computeInvM' in " << dim << " dimension(s)...'" << std::endl;

	//create test problem class object
	//userInputParameters userInputs;
	//userInputs.loadUserInput();

    testInvM<dim,1> test(userInputs);
	//check invM norm
	if ((test.invMNorm - 1700.0) < 1.0e-10) {pass=true;}
	char buffer[100];
	snprintf(buffer, sizeof(buffer), "Test result for 'computeInvM' in   %u dimension(s): %u\n", dim, pass);
	std::cout << buffer;

	return pass;
}
