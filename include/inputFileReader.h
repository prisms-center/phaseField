// Class to read in the user input from parameters.in

#ifndef INCLUDE_INPUTFILEREADER_H_
#define INCLUDE_INPUTFILEREADER_H_

#include <deal.II/base/parameter_handler.h>

#include "variableAttributeLoader.h"

#include <fstream>
#include <string>
#include <vector>

class inputFileReader
{
public:
  // Constructor
  inputFileReader(std::string             input_file_name,
                  variableAttributeLoader variable_attributes);

  // Method to get a list of entry values from multiple subsections in an input
  // file
  std::vector<std::string>
  get_subsection_entry_list(const std::string parameters_file_name,
                            const std::string subsec_name,
                            const std::string entry_name,
                            const std::string default_entry) const;

  // Method to count the number of related entries in an input file
  unsigned int
  get_number_of_entries(const std::string parameters_file_name,
                        const std::string keyword,
                        const std::string entry_name) const;

  // Get the trailing part of the entry name after a specified string (used to
  // extract the model constant names)
  std::vector<std::string>
  get_entry_name_ending_list(const std::string parameters_file_name,
                             const std::string keyword,
                             const std::string entry_name_begining) const;

  // Method to declare the parameters to be read from an input file
  void
  declare_parameters(dealii::ParameterHandler    &parameter_handler,
                     const std::vector<fieldType> var_types,
                     const std::vector<PDEType>   var_eq_types,
                     const unsigned int           num_of_constants,
                     const std::vector<bool>) const;

  // Method to check if a line has the desired contents and if so, extract it
  bool
  parse_line(std::string       line,
             const std::string keyword,
             const std::string entry_name,
             std::string      &out_string,
             bool              expect_equals_sign) const;

  // Variables
  dealii::ParameterHandler parameter_handler;
  std::vector<fieldType>   var_types;
  std::vector<PDEType>     var_eq_types;
  unsigned int             num_pp_vars;
  unsigned int             num_constants;
  std::vector<std::string> model_constant_names;
  std::vector<std::string> var_names;
  unsigned int             number_of_dimensions;
  std::vector<bool>        var_nucleates;
  std::vector<bool>        var_nonlinear;
};

#endif /* INCLUDE_INPUTFILEREADER_H_ */
