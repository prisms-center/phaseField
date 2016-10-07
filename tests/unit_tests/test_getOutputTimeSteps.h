// Unit test(s) for the method "computeInvM"
template <int dim>
class getOutputTimeStepsTest: public MatrixFreePDE<dim>
{
 public: 
	std::vector<int> outputTimeStepList_public;
	getOutputTimeStepsTest(){
    //init the MatrixFreePDE class for testing
    this->initForTests();

    //call computeInvM()
    this->getOutputTimeSteps();

    outputTimeStepList_public = this->outputTimeStepList;
  };

  
 private:
	//RHS implementation for explicit solve
	  void getRHS(const MatrixFree<dim,double> &data,
		      std::vector<vectorType*> &dst,
		      const std::vector<vectorType*> &src,
		      const std::pair<unsigned int,unsigned int> &cell_range) const{};

};

template <int dim,typename T>
  bool unitTest<dim,T>::test_getOutputTimeSteps(){
  	bool pass = true;
	std::cout << "\nTesting 'getOutputTimeSteps' in " << dim << " dimension(s)...'" << std::endl;
 
	//create test problem class object
	getOutputTimeStepsTest<dim> test;

	// Check if calculated value equals expected value
	std::vector<int> expected_result = {0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000};
	for (unsigned int i=0; i < test.outputTimeStepList_public.size(); i++){
		std::cout << test.outputTimeStepList_public[i] << " " << expected_result[i] << std::endl;
		if (test.outputTimeStepList_public[i] != expected_result[i]){
			pass = false;
		}
	}

	char buffer[100];
	sprintf (buffer, "Test result for 'getOutputTimeSteps' in   %u dimension(s): %u\n", dim, pass);
	std::cout << buffer;
	
	return pass;
}
