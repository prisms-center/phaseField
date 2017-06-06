// Class to read in the user input from parameters.in

#ifndef INCLUDE_INPUTFILEREADER_H_
#define INCLUDE_INPUTFILEREADER_H_

#include "dealIIheaders.h"


class inputFileReader
{
public:
	// Method to declare the parameters to be read from an input file
	inputFileReader(dealii::ParameterHandler & parameter_handler, std::string input_file_name, unsigned int number_of_variables, unsigned int number_of_materials, unsigned int number_of_pp_variables);
};

#endif /* INCLUDE_INPUTFILEREADER_H_ */
