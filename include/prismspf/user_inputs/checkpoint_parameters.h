// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>

#include <prismspf/user_inputs/temporal_discretization.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include <climits>
#include <set>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that holds checkpoint parameters.
 */
struct CheckpointParameters
{
  /**
   * @brief Return if the increment should be checkpointted.
   */
  [[nodiscard]] bool
  should_checkpoint(unsigned int increment) const
  {
    return checkpoint_list.contains(increment);
  }

  /**
   * @brief Postprocess and validate parameters.
   */
  void
  validate()
  { // todo
  }

  /**
   * @brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  /**
   * @brief Set the user checkpoint list
   */
  template <typename ListType>
  void
  add_checkpoint_list(const ListType &list)
  {
    checkpoint_list.insert(list.begin(), list.end());
  }

  /**
   * @brief Set the user checkpoint list
   */
  void
  add_equal_spacing_checkpoints(unsigned int num_checkpoints, unsigned int num_increments)
  {
    if (!num_checkpoints)
      {
        return;
      }
    for (unsigned int checkpoint = 0; checkpoint <= num_increments;
         checkpoint += num_increments / num_checkpoints)
      {
        checkpoint_list.insert(checkpoint);
      }
  }

  /**
   * @brief Set the user checkpoint list
   */
  void
  add_log_spacing_checkpoints(unsigned int num_checkpoints, unsigned int num_increments)
  {
    if (!num_checkpoints)
      {
        return;
      }
    for (unsigned int checkpoint = 1; checkpoint <= num_checkpoints; checkpoint++)
      {
        checkpoint_list.insert(static_cast<unsigned int>(
          std::round(std::pow(double(num_increments),
                              double(checkpoint) / double(num_checkpoints)))));
      }
  }

  /**
   * @brief Set the user checkpoint list
   */
  void
  add_n_per_decade_checkpoints(unsigned int num_checkpoints, unsigned int num_increments)
  {
    if (!num_checkpoints)
      {
        return;
      }
    AssertThrow(num_increments > 1,
                dealii::ExcMessage("For n per decaded spaced checkpoints, the number of "
                                   "increments must be greater than 1."));

    for (unsigned int iteration = 2; iteration <= num_increments; iteration++)
      {
        const auto decade = static_cast<unsigned int>(std::ceil(std::log10(iteration)));
        const auto step_size =
          static_cast<unsigned int>(std::pow(10, decade) / num_checkpoints);
        if (iteration % step_size == 0)
          {
            checkpoint_list.insert(iteration);
          }
      }
  }

  /**
   * @brief Set the user checkpoint list
   */
  void
  clear_checkpoint_list()
  {
    checkpoint_list.clear();
  }

  /**
   * @brief Get the number of checkpoints that will be made
   */
  [[nodiscard]] unsigned int
  get_num_checkpoints() const
  {
    return checkpoint_list.size();
  }

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler) const;

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  void
  assign_parameters(dealii::ParameterHandler     &parameter_handler,
                    const TemporalDiscretization &temporal_discretization);

  // Whether to load from a checkpoint
  bool should_load_checkpoint = false;

  // Checkpoint file name
  std::string file_name;

  // Whether to print timing information with checkpoint
  // TODO (landinjm): Implement this.
  bool print_timing_with_checkpoint = false;

  // List of increments that checkpoint the solution to file
  std::set<unsigned int> checkpoint_list;
};

inline void
CheckpointParameters::print_parameter_summary() const
{
  ConditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Checkpoint Parameters\n"
    << "================================================\n"
    << "Checkpoint file name: " << file_name << "\n"
    << "Number of checkpoints: " << get_num_checkpoints() << "\n"
    << "Print timing info: " << bool_to_string(print_timing_with_checkpoint) << "\n";

  ConditionalOStreams::pout_summary() << "Checkpoint increment list: ";
  for (const auto &iteration : checkpoint_list)
    {
      ConditionalOStreams::pout_summary() << iteration << " ";
    }
  ConditionalOStreams::pout_summary() << "\n\n" << std::flush;
}

inline void
CheckpointParameters::declare_parameters(
  dealii::ParameterHandler &parameter_handler) const
{
  parameter_handler.enter_subsection("checkpoints");
  {
    parameter_handler.declare_entry(
      "load from checkpoint",
      "false",
      dealii::Patterns::Bool(),
      "Whether to load from a checkpoint created during a previous simulation.");
    parameter_handler.declare_entry(
      "condition",
      "EQUAL_SPACING",
      dealii::Patterns::Selection("EQUAL_SPACING|LOG_SPACING|N_PER_DECADE|LIST"),
      "The spacing type for saving checkpoints (either EQUAL_SPACING, "
      "LOG_SPACING, N_PER_DECADE, or LIST).");
    parameter_handler.declare_entry(
      "list",
      "0",
      dealii::Patterns::List(dealii::Patterns::Integer(0, INT_MAX), 0, INT_MAX, ","),
      "The list of time steps to save checkpoints, used for the LIST type.");
    parameter_handler.declare_entry(
      "number",
      "0",
      dealii::Patterns::Integer(0, INT_MAX),
      "The number of checkpoints (or number of checkpoints per decade for the "
      "N_PER_DECADE type).");
  }
  parameter_handler.leave_subsection();
}

inline void
CheckpointParameters::assign_parameters(
  dealii::ParameterHandler     &parameter_handler,
  const TemporalDiscretization &temporal_discretization)
{
  parameter_handler.enter_subsection("checkpoints");
  {
    should_load_checkpoint = parameter_handler.get_bool("load from checkpoint");

    std::string  condition = parameter_handler.get("condition");
    unsigned int num_checkpoints =
      static_cast<unsigned int>(parameter_handler.get_integer("number"));
    unsigned int num_increments = temporal_discretization.num_increments;
    if (condition == "EQUAL_SPACING")
      {
        add_equal_spacing_checkpoints(num_checkpoints, num_increments);
      }
    else if (condition == "LOG_SPACING")
      {
        add_log_spacing_checkpoints(num_checkpoints, num_increments);
      }
    else if (condition == "N_PER_DECADE")
      {
        add_n_per_decade_checkpoints(num_checkpoints, num_increments);
      }
    add_checkpoint_list(dealii::Utilities::string_to_int(
      dealii::Utilities::split_string_list(parameter_handler.get("list"))));
  }
  parameter_handler.leave_subsection();
}

PRISMS_PF_END_NAMESPACE
