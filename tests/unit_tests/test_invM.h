// Unit test(s) for the method "computeInvM"
template <int dim, int degree>
class testInvM: public MatrixFreePDE<dim,degree>
{
 public:
  testInvM(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) {

	  // Initialize the test field object (needed for computeInvM())
	  Field<dim> test_field(SCALAR,PARABOLIC,"c");
	  this->fields.push_back(test_field);

	  //init the MatrixFreePDE class for testing
	  this->initForTests();

	  //call computeInvM()
	  this->computeInvM();
	  invMNorm=this->invM.l2_norm();

	  // Need to clear fields or there's an error in the destructor
	  this->fields.clear();
  };
  double invMNorm;

  void setBCs(){};

 private:
  //RHS implementation for explicit solve
  void getRHS(const MatrixFree<dim,double> &data,
	      std::vector<vectorType*> &dst,
	      const std::vector<vectorType*> &src,
	      const std::pair<unsigned int,unsigned int> &cell_range) const{};

  void residualRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
  					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {};

  void residualLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
  					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {};

  void energyDensity(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list, const dealii::VectorizedArray<double> & JxW_value,
  			  	  	 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {};

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
	sprintf (buffer, "Test result for 'computeInvM' in   %u dimension(s): %u\n", dim, pass);
	std::cout << buffer;

	return pass;
}
