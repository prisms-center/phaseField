// get_output_timesteps
// This function generates a list of time steps where the program should generate output. The inputs
// are read in from the file "parameters.h" via #define macros.

template <int dim>
void MatrixFreePDE<dim>::getOutputTimeSteps(){

	// Determine the maximum number of increments
	unsigned int maxIncrements;
	if (std::ceil(timeFinal/timeStep) < timeIncrements){
		maxIncrements = std::ceil(timeFinal/timeStep);
	}
	else {
		maxIncrements = timeIncrements;
	}

	if (outputCondition == "EQUAL_SPACING"){
		for (unsigned int iter = 0; iter <= maxIncrements; iter += maxIncrements/numOutputs){
			outputTimeStepList.push_back(iter);
			std::cout << "output iter: " << iter << std::endl;
		}
	}
	else if (outputCondition == "LOG_SPACING"){
		outputTimeStepList.push_back(0);
		for (unsigned int output = 1; output <= numOutputs; output++){
			outputTimeStepList.push_back(std::round(std::pow(10,double(output)/double(numOutputs)*std::log10(maxIncrements))));
			std::cout << "output iter: " << std::round(std::pow(10,double(output)/double(numOutputs)*std::log10(maxIncrements))) << std::endl;
		}
	}
	else if (outputCondition == "N_PER_DECADE"){
		outputTimeStepList.push_back(0);
		outputTimeStepList.push_back(1);
		for (unsigned int iter = 2; iter <= maxIncrements; iter++){
			int decade = std::ceil(std::log10(iter));
			int step_size = (std::pow(10,decade))/numOutputs;
			if (iter%step_size == 0){
				outputTimeStepList.push_back(iter);
				std::cout << "output iter: " << iter << std::endl;
			}

		}
	}




}




