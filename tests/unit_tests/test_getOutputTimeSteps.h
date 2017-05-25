// Unit test(s) for the method "getOutputTimeSteps"
template <int dim, int degree>
class getOutputTimeStepsTest: public MatrixFreePDE<dim,degree>
{
	public:
	std::vector<unsigned int> outputTimeStepList_public;
	getOutputTimeStepsTest(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs){
	};


	void getTimeStepList(std::string outputSpacingType, unsigned int numberOfOutputs, std::vector<unsigned int> userGivenTimeStepList){
		this->getOutputTimeSteps(outputSpacingType,numberOfOutputs,userGivenTimeStepList,this->outputTimeStepList);
		outputTimeStepList_public = this->outputTimeStepList;
	};


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
bool unitTest<dim,T>::test_getOutputTimeSteps(std::string outputSpacingType, unsigned int numberOfOutputs, std::vector<unsigned int> userGivenTimeStepList, userInputParameters<dim> userInputs){
	bool pass = true;
	std::cout << "\nTesting 'getOutputTimeSteps' for type " << outputSpacingType << " ...'" << std::endl;

	//create test problem class object
	getOutputTimeStepsTest<dim,finiteElementDegree> test(userInputs);

	test.getTimeStepList(outputSpacingType,numberOfOutputs,userGivenTimeStepList);

	// Check if calculated value equals expected value
	std::vector<unsigned int> expected_result;
	if (outputSpacingType == "EQUAL_SPACING"){
		expected_result = {0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000};
	}
	else if (outputSpacingType == "LOG_SPACING"){
		expected_result = {0,3,7,20,53,141,381,1025,2759,7429,20000};
	}
	else if (outputSpacingType == "N_PER_DECADE"){
		expected_result = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
				200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000,
				8000, 9000, 10000, 20000};
	}
	else if (outputSpacingType == "LIST"){
			expected_result = {0, 3, 55, 61};
	}

	for (unsigned int i=0; i < test.outputTimeStepList_public.size(); i++){
		if (test.outputTimeStepList_public[i] != expected_result[i]){
			pass = false;
		}
	}

	char buffer[100];
	//char temp[] = outputSpacingType.c_str();
	sprintf (buffer, "Test result for 'getOutputTimeSteps' for type %s : %u\n", outputSpacingType.c_str(), pass);
	std::cout << buffer;
	
	return pass;
}
