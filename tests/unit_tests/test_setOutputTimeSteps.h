// Unit test(s) for the method "setOutputTimeSteps"

template <int dim, typename T>
bool
unitTest<dim, T>::test_setOutputTimeSteps()
{
  bool pass = false;
  char buffer[100];
  std::cout << "\nTesting 'setOutputTimeSteps' via the public "
               "'loadInputParameters' method...'"
            << std::endl;

  dealii::ParameterHandler parameter_handler;
  variableAttributeLoader  attribute_loader;
  const AttributesList     var_attributes = attribute_loader.get_var_attributes();
  const AttributesList     pp_attributes  = attribute_loader.get_pp_attributes();
  inputFileReader input_file_reader("parameters_test.prm", var_attributes, pp_attributes);
  input_file_reader.declare_parameters(parameter_handler, 0);
#if (DEAL_II_VERSION_MAJOR < 9 && DEAL_II_VERSION_MINOR < 5)
  parameter_handler.read_input("parameters_test.prm");
#else
  parameter_handler.parse_input("parameters_test.prm");
#endif

  // userInputParameters<dim> userInputs;

  // Subtest 1 (EQUAL_SPACING)
  parameter_handler.set("Output condition", "EQUAL_SPACING");

  // userInputs.loadInputParameters(input_file_reader,parameter_handler,var_types.size(),0,0,0);
  userInputParameters<dim> userInputs1(input_file_reader,
                                       parameter_handler,
                                       var_attributes,
                                       pp_attributes);

  std::vector<unsigned int> expected_result =
    {0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000};

  bool pass1 = true;
  for (unsigned int i = 0; i < userInputs1.outputTimeStepList.size(); i++)
    {
      if (userInputs1.outputTimeStepList[i] != expected_result[i])
        {
          pass1 = false;
        }
    }

  snprintf(buffer,
           sizeof(buffer),
           "Subtest 1 for 'setOutputTimeSteps' (EQUAL_SPACING) result: %u\n",
           pass1);
  std::cout << buffer;

  // Subtest 2 (LOG_SPACING)
  parameter_handler.set("Output condition", "LOG_SPACING");
  userInputParameters<dim> userInputs2(input_file_reader,
                                       parameter_handler,
                                       var_attributes,
                                       pp_attributes);

  expected_result = {0, 3, 7, 20, 53, 141, 381, 1025, 2759, 7429, 20000};

  bool pass2 = true;
  for (unsigned int i = 0; i < userInputs2.outputTimeStepList.size(); i++)
    {
      if (userInputs2.outputTimeStepList[i] != expected_result[i])
        {
          pass2 = false;
        }
    }

  snprintf(buffer,
           sizeof(buffer),
           "Subtest 2 for 'setOutputTimeSteps' (LOG_SPACING) result: %u\n",
           pass2);
  std::cout << buffer;

  // Subtest 3 (N_PER_DECADE)
  parameter_handler.set("Output condition", "N_PER_DECADE");
  userInputParameters<dim> userInputs3(input_file_reader,
                                       parameter_handler,
                                       var_attributes,
                                       pp_attributes);

  expected_result = {0,    1,    2,    3,    4,    5,    6,    7,     8,    9,
                     10,   20,   30,   40,   50,   60,   70,   80,    90,   100,
                     200,  300,  400,  500,  600,  700,  800,  900,   1000, 2000,
                     3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000};

  bool pass3 = true;
  for (unsigned int i = 0; i < userInputs3.outputTimeStepList.size(); i++)
    {
      if (userInputs3.outputTimeStepList[i] != expected_result[i])
        {
          pass3 = false;
        }
    }

  snprintf(buffer,
           sizeof(buffer),
           "Subtest 3 for 'setOutputTimeSteps' (N_PER_DECADE) result: %u\n",
           pass3);
  std::cout << buffer;

  // Subtest 4 (LIST)
  parameter_handler.set("Output condition", "LIST");
  userInputParameters<dim> userInputs4(input_file_reader,
                                       parameter_handler,
                                       var_attributes,
                                       pp_attributes);

  expected_result = {0, 3, 55, 61};

  bool pass4 = true;
  for (unsigned int i = 0; i < userInputs4.outputTimeStepList.size(); i++)
    {
      if (userInputs4.outputTimeStepList[i] != expected_result[i])
        {
          pass4 = false;
        }
    }

  snprintf(buffer,
           sizeof(buffer),
           "Subtest 4 for 'setOutputTimeSteps' (LIST) result: %u\n",
           pass4);
  std::cout << buffer;

  if (pass4 && pass3 && pass2 && pass1)
    {
      pass = true;
    }

  snprintf(buffer, sizeof(buffer), "Test result for 'setOutputTimeSteps': %u\n", pass);
  std::cout << buffer;

  return pass;
}
