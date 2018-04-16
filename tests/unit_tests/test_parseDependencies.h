template <int dim>
class testParseDependencies: public userInputParameters<dim>
{
 public:
     // This is supposed to be a stripped down test class for userInputParameters
     // Note: this doesn't work
     testParseDependencies(inputFileReader, dealii::ParameterHandler, variableAttributeLoader): userInputParameters<dim>(inputFileReader, dealii::ParameterHandler, variableAttributeLoader) {};

     void tester(){
         result = "this is a test!";
     };

     std::string result;
};



template <int dim,typename T>
bool unitTest<dim,T>::test_parseDependencies(){
    char buffer[100];

	std::cout << "\nTesting 'parseDependencies'... " << std::endl;

    variableAttributeLoader variable_attributes;
    inputFileReader input_file_reader("parameters_test.in",variable_attributes);
    userInputParameters<2> userInputs(input_file_reader,input_file_reader.parameter_handler,variable_attributes);

    testParseDependencies<2> test(input_file_reader,input_file_reader.parameter_handler,variable_attributes);
    test.tester();

    std::cout << test.result << std::endl;

}
