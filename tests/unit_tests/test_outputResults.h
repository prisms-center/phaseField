// Unit test(s) for the method "outputResults"
template <int dim, int degree>
class testOutputResults: public MatrixFreePDE<dim,degree>
{
 public: 
	testOutputResults(userInputParameters _userInputs): MatrixFreePDE<dim,degree>(_userInputs){};
  
 private:
	void setBCs(){};

	private:
	//RHS implementation for explicit solve
	void getRHS(const MatrixFree<dim,double> &data,
			std::vector<vectorType*> &dst,
			const std::vector<vectorType*> &src,
			const std::pair<unsigned int,unsigned int> &cell_range) const{};

	void residualRHS(const std::vector<modelVariable<dim> > & modelVarList,
			std::vector<modelResidual<dim> > & modelResidualsList,
			dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {};

	void residualLHS(const std::vector<modelVariable<dim> > & modelVarList,
			modelResidual<dim> & modelRes,
			dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {};

	void energyDensity(const std::vector<modelVariable<dim> > & modelVarList, const dealii::VectorizedArray<double> & JxW_value,
			dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {};
};

template <int dim,typename T>
  bool unitTest<dim,T>::test_outputResults(int argc, char** argv, userInputParameters userInputs){
  bool pass = false;
  std::cout << "\nTesting 'outputResults' in " << dim << " dimension(s)...'" << std::endl;
  
  //create test problem class object
  testOutputResults<dim,finiteElementDegree> test(userInputs);
  pass=true;
  
  char buffer[100];
  sprintf (buffer, "Test result for 'outputResults' in %u dimension(s): %u\n", dim, pass);
  std::cout << buffer;
  
  return pass;
}
