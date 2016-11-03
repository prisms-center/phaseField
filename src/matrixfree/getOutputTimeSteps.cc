// get_output_timesteps
// This function generates a list of time steps where the program should generate output. The inputs
// are read in from the file "parameters.h" via #define macros.

template <int dim>
void MatrixFreePDE<dim>::getOutputTimeSteps(std::string outputSpacingType, unsigned int numberOfOutputs, std::vector<unsigned int> & userGivenTimeStepList, std::vector<unsigned int> & timeStepList){

	// Determine the maximum number of increments
	if (outputSpacingType == "EQUAL_SPACING"){
		if (numberOfOutputs > totalIncrements)
			numberOfOutputs = totalIncrements;

		for (unsigned int iter = 0; iter <= totalIncrements; iter += totalIncrements/numberOfOutputs){
			timeStepList.push_back(iter);
		}
	}
	else if (outputSpacingType == "LOG_SPACING"){
		timeStepList.push_back(0);
		for (unsigned int output = 1; output <= numberOfOutputs; output++){
			timeStepList.push_back(std::round(std::pow(10,double(output)/double(numberOfOutputs)*std::log10(totalIncrements))));		}
	}
	else if (outputSpacingType == "N_PER_DECADE"){
		timeStepList.push_back(0);
		timeStepList.push_back(1);
		for (unsigned int iter = 2; iter <= totalIncrements; iter++){
			int decade = std::ceil(std::log10(iter));
			int step_size = (std::pow(10,decade))/numberOfOutputs;
			if (iter%step_size == 0){
				timeStepList.push_back(iter);
			}

		}
	}
	else if (outputSpacingType == "LIST"){
		timeStepList = userGivenTimeStepList;
	}
}




