// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/config.h>

#include <map>
#include <set>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

struct variableAttributes;

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
  inputFileReader(std::string                                       input_file_name,
                  const std::map<unsigned int, variableAttributes> &_var_attributes);

  /**
   * \brief Get the trailing part of the entry name after a specified string (used to
   * extract the model constant names).
   */
  [[nodiscard]] std::set<std::string>
  get_model_constant_names();

  /**
   * \brief Get the variable attributes.
   */
  [[nodiscard]] std::map<unsigned int, variableAttributes>
  get_var_attributes() const
  {
    return var_attributes;
  }

  /**
   * \brief Get the parameter handler.
   */
  [[nodiscard]] dealii::ParameterHandler &
  get_parameter_handler()
  {
    return parameter_handler;
  }

  /**
   * \brief Get the model constant names.
   */
  [[nodiscard]] const std::set<std::string> &
  get_model_constant_names() const
  {
    return model_constant_names;
  }

  /**
   * \brief Get the number of dimensions.
   */
  [[nodiscard]] unsigned int
  get_dim() const
  {
    return number_of_dimensions;
  };

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
  check_keyword_match(const std::string &line, const std::string &keyword);

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

private:
  std::string                                       parameters_file_name;
  const std::map<unsigned int, variableAttributes> &var_attributes;
  dealii::ParameterHandler                          parameter_handler;
  std::set<std::string>                             model_constant_names;
  unsigned int                                      number_of_dimensions = 0;
};

PRISMS_PF_END_NAMESPACE
