#ifndef input_file_reader_h
#define input_file_reader_h

#include <deal.II/base/parameter_handler.h>

#include <core/variable_attribute_loader.h>
#include <set>
#include <string>
#include <vector>

/**
 * \brief Parameters file reader. Declares parameter names in a dealii parameter_handler
 * and parses the file for the values. Variable assignment occurs in userInputParameters.
 */
class inputFileReader
{
public:
  /**
   * \brief Constructor.
   */
  inputFileReader(const std::string    &input_file_name,
                  const AttributesList &_var_attributes,
                  const AttributesList &_pp_attributes);

  /**
   * \brief Method to get a list of entry values from multiple subsections in an input
   * file.
   */
  [[nodiscard]] std::vector<std::string>
  get_subsection_entry_list(const std::string &subsec_name,
                            const std::string &entry_name,
                            const std::string &default_entry);

  /**
   * \brief Get the trailing part of the entry name after a specified string (used to
   * extract the model constant names).
   */
  [[nodiscard]] std::set<std::string>
  get_model_constant_names();

  /**
   * \brief Method to declare the parameters to be read from an input file.
   */
  void
  declare_parameters();

  /**
   * \brief Method to check if a line has the desired contents and if so, extract it.
   */
  bool
  parse_line(std::string        line,
             const std::string &keyword,
             const std::string &entry_name,
             std::string       &out_string,
             bool               expect_equals_sign);

  /**
   * \brief Strip spaces from the front and back of a string.
   */
  void
  strip_spaces(std::string &line);

  /**
   * \brief Check whether a string starts with a keyword.
   */
  bool
  check_keyword_match(std::string &line, const std::string &keyword);

  /**
   * \brief Declare parameters for the mesh.
   */
  void
  declare_mesh();

  /**
   * \brief Declare parameters for timestepping.
   */
  void
  declare_time_discretization();

  /**
   * \brief Declare parameters for linear and nonlinear solvers.
   */
  void
  declare_solver_parameters();

  /**
   * \brief Declare parameters for outputs.
   */
  void
  declare_output_parameters();

  /**
   * \brief Declare parameters for loading ICs from files.
   */
  void
  declare_load_IC_parameters();

  /**
   * \brief Declare parameters for checkpoints.
   */
  void
  declare_checkpoint_parameters();

  /**
   * \brief Declare parameters for boundary conditions.
   */
  void
  declare_BC_parameters();

  /**
   * \brief Declare parameters for pinned points.
   */
  void
  declare_pinning_parameters();

  /**
   * \brief Declare parameters for nucleation
   */
  void
  declare_nucleation_parameters();

  /**
   * \brief Declare parameters for grain remapping.
   */
  void
  declare_grain_remapping_parameters();

  /**
   * \brief Declare parameters for grain structure loading.
   */
  void
  declare_grain_loading_parameters();

  /**
   * \brief Declare parameters for user-defined model constants.
   */
  void
  declare_model_constants();

  const std::string        parameters_file_name;
  const AttributesList    &var_attributes;
  const AttributesList    &pp_attributes;
  dealii::ParameterHandler parameter_handler;
  std::set<std::string>    model_constant_names;
  uint                     number_of_dimensions;
};

#endif /* INCLUDE_INPUTFILEREADER_H_ */
