// Class to read in the user input from parameters.in

#ifndef INCLUDE_INPUTFILEREADER_H_
#define INCLUDE_INPUTFILEREADER_H_

#include "dealIIheaders.h"
#include <fstream>

class inputFileReader
{
public:
	// Method to get a list of entry values from multiple subsections in an input file
	std::vector<std::string> get_subsection_entry_list(const std::string parameters_file_name, const std::string subsec_name,
														const std::string entry_name, const std::string default_entry) const;

	// Method to count the number of related entries in an input file
	unsigned int get_number_of_entries(const std::string parameters_file_name,const std::string keyword, const std::string entry_name) const;

	// Method to extract the value of a single parameter from a deal.II-style input file
	std::string get_last_value_of_parameter(const std::string &parameters, const std::string &parameter_name) const;

	// Method to declare the parameters to be read from an input file
	void declare_parameters(dealii::ParameterHandler & parameter_handler,
							const std::vector<std::string> var_types, const unsigned int number_of_materials, const unsigned int number_of_pp_variables,
							const unsigned int num_of_constants) const;

	// Method to check if a line has the desired contents and if so, extract it
	bool parse_line(std::string line, const std::string keyword, const std:: string entry_name, std::string & out_string, bool expect_equals_sign) const;

};

#endif /* INCLUDE_INPUTFILEREADER_H_ */
