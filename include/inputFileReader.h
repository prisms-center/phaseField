// Class to read in the user input from parameters.in

#ifndef INCLUDE_INPUTFILEREADER_H_
#define INCLUDE_INPUTFILEREADER_H_

#include <deal.II/base/parameter_handler.h>

#include "variableAttributes.h"

#include <string>
#include <vector>

class inputFileReader
{
public:
  // Constructor
  inputFileReader(const std::string    &input_file_name,
                  const AttributesList &_var_attributes,
                  const AttributesList &_pp_attributes);

  // Method to get a list of entry values from multiple subsections in an input
  // file
  std::vector<std::string>
  get_subsection_entry_list(const std::string &parameters_file_name,
                            const std::string &subsec_name,
                            const std::string &entry_name,
                            const std::string &default_entry) const;

  // Method to count the number of related entries in an input file
  unsigned int
  get_number_of_entries(const std::string &parameters_file_name,
                        const std::string &keyword,
                        const std::string &entry_name) const;

  // Get the trailing part of the entry name after a specified string (used to
  // extract the model constant names)
  std::vector<std::string>
  get_entry_name_ending_list(const std::string &parameters_file_name,
                             const std::string &keyword,
                             const std::string &entry_name_begining) const;

  // Method to declare the parameters to be read from an input file
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler,
                     const unsigned int        num_of_constants) const;

  // Method to check if a line has the desired contents and if so, extract it
  bool
  parse_line(std::string        line,
             const std::string &keyword,
             const std::string &entry_name,
             std::string       &out_string,
             bool               expect_equals_sign) const;

  // Variables
  const AttributesList    &var_attributes;
  const AttributesList    &pp_attributes;
  dealii::ParameterHandler parameter_handler;
  unsigned int             num_constants;
  std::vector<std::string> model_constant_names;
  unsigned int             number_of_dimensions;
};

#endif /* INCLUDE_INPUTFILEREADER_H_ */
