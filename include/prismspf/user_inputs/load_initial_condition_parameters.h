// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/types.h>

#include <boost/algorithm/string/predicate.hpp>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that store the read-in information for a single file
 */
struct InitialConditionFile
{
  // File name
  std::string filename;

  // Grid type
  DataFormatType dataset_format;

  // File variable names
  std::vector<std::string> file_variable_names;

  // Simulation variable names
  std::vector<std::string> simulation_variable_names;

  // Number of data points in each direction
  std::array<dealii::types::global_dof_index, 3> n_data_points = {
    {0, 0, 0}
  };
};

/**
 * @brief Struct that stores relevant load initial condition information
 */
struct LoadInitialConditionParameters
{
public:
  /**
   * @brief Maximum number of initial condition files
   */
  static constexpr unsigned int max_files = 8;

  /**
   * @brief Postprocess and validate parameters.
   */
  void
  postprocess_and_validate();

  /**
   * @brief Clear the initial condition parameters.
   */
  void
  clear()
  {
    read_initial_conditions_from_file = false;
    ic_files.clear();
  }

  /**
   * @brief Set the read initial conditions from file flag.
   */
  void
  set_read_initial_conditions_from_file(bool _read_initial_conditions_from_file)
  {
    read_initial_conditions_from_file = _read_initial_conditions_from_file;
  }

  /**
   * @brief Get the read initial conditions from file flag.
   */
  [[nodiscard]] bool
  get_read_initial_conditions_from_file() const
  {
    return read_initial_conditions_from_file;
  }

  /**
   * @brief Add a initial condition file.
   */
  void
  add_initial_condition_file(InitialConditionFile _ic_file)
  {
    if (read_initial_conditions_from_file)
      {
        ic_files.push_back(_ic_file);
      }
  }

  /**
   * @brief Get the number of initial condition files.
   */
  [[nodiscard]] unsigned int
  get_n_initial_condition_files() const
  {
    return ic_files.size();
  }

  /**
   * @brief Get the initial condition files.
   */
  [[nodiscard]] const std::vector<InitialConditionFile> &
  get_initial_condition_files() const
  {
    return ic_files;
  }

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler,
                     unsigned int              max_criteria = 5) const;

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  void
  assign_parameters(dealii::ParameterHandler &parameter_handler,
                    unsigned int              max_criteria = 5);

private:
  // Whether to read initial conditions from file
  bool read_initial_conditions_from_file = false;

  // IC files
  std::vector<InitialConditionFile> ic_files;
};

inline void
LoadInitialConditionParameters::postprocess_and_validate()
{
  for (const auto &ic_file : ic_files)
    {
      // Check that the file variables are the same length as the simulation variables
      AssertThrow(ic_file.file_variable_names.size() ==
                    ic_file.simulation_variable_names.size(),
                  dealii::ExcMessage("The number of file variables must be the same as "
                                     "the number of simulation variables"));
    }

  // TODO (landinjm): Check that there are no duplicate field names so we don't double
  // assign.
}

inline void
LoadInitialConditionParameters::declare_parameters(
  dealii::ParameterHandler &parameter_handler,
  unsigned int              max_criteria) const
{
  parameter_handler.declare_entry("read initial conditions from file",
                                  "false",
                                  dealii::Patterns::Bool(),
                                  "Whether to read any initial conditions from file.");

  for (unsigned int i = 0; i < max_criteria; i++)
    {
      parameter_handler.enter_subsection("initial condition file " + std::to_string(i));
      {
        parameter_handler.declare_entry("file name",
                                        "",
                                        dealii::Patterns::Anything(),
                                        "The file name to load from for each variable.");
        parameter_handler.declare_entry("dataset format",
                                        "vtk_unstructured_grid",
                                        dealii::Patterns::Anything(),
                                        "The type of grid in the file.");
        parameter_handler.declare_entry("file variable names",
                                        "",
                                        dealii::Patterns::List(
                                          dealii::Patterns::Anything()),
                                        "The name of the variable in the file.");
        parameter_handler.declare_entry("simulation variable names",
                                        "",
                                        dealii::Patterns::List(
                                          dealii::Patterns::Anything()),
                                        "The name of the variable in the file.");
        parameter_handler.declare_entry(
          "data points in x direction",
          "-1",
          dealii::Patterns::Integer(-1, INT_MAX),
          "The number of data points of the input file in the x direction.");
        parameter_handler.declare_entry(
          "data points in y direction",
          "-1",
          dealii::Patterns::Integer(-1, INT_MAX),
          "The number of data points of the input file in the y direction.");
        parameter_handler.declare_entry(
          "data points in z direction",
          "-1",
          dealii::Patterns::Integer(-1, INT_MAX),
          "The number of data points of the input file in the z direction.");
      }
      parameter_handler.leave_subsection();
    }
}

inline void
LoadInitialConditionParameters::assign_parameters(
  dealii::ParameterHandler &parameter_handler,
  unsigned int              max_criteria)
{
  set_read_initial_conditions_from_file(
    parameter_handler.get_bool("read initial conditions from file"));
  static std::array<std::string, 3> axis_labels = {
    {"x", "y", "z"}
  };
  for (unsigned int i = 0; i < max_criteria; i++)
    {
      parameter_handler.enter_subsection("initial condition file " + std::to_string(i));
      {
        // Check if the file is specified
        if (!parameter_handler.get("file name").empty())
          {
            // Create the LoadICFile object
            InitialConditionFile ic_file;
            ic_file.filename              = parameter_handler.get("file name");
            const std::string type_string = parameter_handler.get("dataset format");
            bool              found_type  = false;
            for (unsigned int j = 0;
                 j < static_cast<unsigned int>(DataFormatType::LastEntry);
                 j++)
              {
                if (boost::iequals(type_string,
                                   to_string(static_cast<DataFormatType>(j))))
                  {
                    ic_file.dataset_format = static_cast<DataFormatType>(j);
                    found_type             = true;
                    break;
                  }
              }
            AssertThrow(found_type,
                        dealii::ExcMessage("Unsupported dataset format: " + type_string));
            ic_file.file_variable_names = dealii::Utilities::split_string_list(
              parameter_handler.get("file variable names"));
            ic_file.simulation_variable_names = dealii::Utilities::split_string_list(
              parameter_handler.get("simulation variable names"));
            // Defaults to 0 for unused dimensions/cases that don't require it
            for (unsigned int k = 0; k < 3; ++k)
              {
                ic_file.n_data_points.at(k) =
                  static_cast<unsigned int>(parameter_handler.get_integer(
                    "data points in " + axis_labels.at(k) + " direction"));
              }
            add_initial_condition_file(ic_file);
          }
      }
      parameter_handler.leave_subsection();
    }
}

PRISMS_PF_END_NAMESPACE
