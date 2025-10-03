// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/types.h>

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
   * @brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

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
LoadInitialConditionParameters::print_parameter_summary() const
{
  if (read_initial_conditions_from_file)
    {
      ConditionalOStreams::pout_summary()
        << "================================================\n"
        << "  Load IC Parameters\n"
        << "================================================\n";

      for (const auto &ic_file : ic_files)
        {
          ConditionalOStreams::pout_summary()
            << "File name: " << ic_file.filename << "\n"
            << "Dataset format: " << to_string(ic_file.dataset_format) << "\n"
            << "File variable names: ";
          for (const auto &file_variable_name : ic_file.file_variable_names)
            {
              ConditionalOStreams::pout_summary() << file_variable_name << " ";
            }
          ConditionalOStreams::pout_summary() << "\n"
                                              << "Simulation variable names: ";
          for (const auto &simulation_variable_name : ic_file.simulation_variable_names)
            {
              ConditionalOStreams::pout_summary() << simulation_variable_name << " ";
            }
          if (ic_file.dataset_format == FlatBinary)
            {
              ConditionalOStreams::pout_summary() << "\n Data points in each direction: ";
              for (const auto &n_data_points : ic_file.n_data_points)
                {
                  ConditionalOStreams::pout_summary() << n_data_points << " ";
                }
            }
        }

      ConditionalOStreams::pout_summary() << "\n" << std::flush;
    }
}

PRISMS_PF_END_NAMESPACE