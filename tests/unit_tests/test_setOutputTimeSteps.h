// Unit test(s) for the method "setOutputTimeSteps"

template <int dim,typename T>
bool unitTest<dim,T>::test_setOutputTimeSteps(){
	bool pass = true;
    char buffer[100];
	std::cout << "\nTesting 'setOutputTimeSteps' via the public 'loadInputParameters' method...'" << std::endl;

    dealii::ParameterHandler parameter_handler;
    inputFileReader input_file_reader;
    std::vector<std::string> var_types;
    var_types.push_back("SCALAR");
    input_file_reader.declare_parameters(parameter_handler,var_types,0,0,0);
    parameter_handler.read_input("parameters_test.in");

    userInputParameters<dim> userInputs;

    // Subtest 1 (EQUAL_SPACING)
    parameter_handler.set("Output condition","EQUAL_SPACING");
    userInputs.loadInputParameters(parameter_handler,var_types.size(),0,0,0);

    std::vector<unsigned int> expected_result = {0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000};

    bool pass1 = true;
    for (unsigned int i=0; i < userInputs.outputTimeStepList.size(); i++){
		if (userInputs.outputTimeStepList[i] != expected_result[i]){
			pass1 = false;
		}
	}

    sprintf (buffer, "Subtest 1 for 'setOutputTimeSteps' (EQUAL_SPACING) result: %u\n", pass1);
	std::cout << buffer;

    // Subtest 2 (LOG_SPACING)
    parameter_handler.set("Output condition","LOG_SPACING");
    userInputs.loadInputParameters(parameter_handler,var_types.size(),0,0,0);

    expected_result = {0,3,7,20,53,141,381,1025,2759,7429,20000};

    bool pass2 = true;
    for (unsigned int i=0; i < userInputs.outputTimeStepList.size(); i++){
		if (userInputs.outputTimeStepList[i] != expected_result[i]){
			pass2 = false;
		}
	}

    sprintf (buffer, "Subtest 2 for 'setOutputTimeSteps' (LOG_SPACING) result: %u\n", pass2);
	std::cout << buffer;

    // Subtest 3 (N_PER_DECADE)
    parameter_handler.set("Output condition","N_PER_DECADE");
    userInputs.loadInputParameters(parameter_handler,var_types.size(),0,0,0);

    expected_result = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
            200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000,
            8000, 9000, 10000, 20000};

    bool pass3 = true;
    for (unsigned int i=0; i < userInputs.outputTimeStepList.size(); i++){
		if (userInputs.outputTimeStepList[i] != expected_result[i]){
			pass3 = false;
		}
	}

    sprintf (buffer, "Subtest 3 for 'setOutputTimeSteps' (N_PER_DECADE) result: %u\n", pass3);
	std::cout << buffer;

    // Subtest 4 (LIST)
    parameter_handler.set("Output condition","LIST");
    userInputs.loadInputParameters(parameter_handler,var_types.size(),0,0,0);

    expected_result = {0, 3, 55, 61};

    bool pass4 = true;
    for (unsigned int i=0; i < userInputs.outputTimeStepList.size(); i++){
		if (userInputs.outputTimeStepList[i] != expected_result[i]){
			pass4 = false;
		}
	}

    sprintf (buffer, "Subtest 4 for 'setOutputTimeSteps' (LIST) result: %u\n", pass4);
	std::cout << buffer;


    if (pass4&&pass3&&pass2&&pass1) {pass=true;}

	sprintf (buffer, "Test result for 'setOutputTimeSteps': %u\n", pass);
	std::cout << buffer;

	return pass;
}
