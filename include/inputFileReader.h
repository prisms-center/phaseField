// Class to read in the user input from parameters.in

#ifndef INCLUDE_INPUTFILEREADER_H_
#define INCLUDE_INPUTFILEREADER_H_

#include "dealIIheaders.h"


class inputFileReader
{
public:
	// Method to declare the parameters to be read from an input file
	inputFileReader(dealii::ParameterHandler & parameter_handler, std::string input_file_name);
};

#endif /* INCLUDE_INPUTFILEREADER_H_ */
