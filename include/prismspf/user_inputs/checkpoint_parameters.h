// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

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
public:
  /**
   * @brief Set whether to load from a checkpoint file.
   */
  void
  set_should_load_checkpoint(bool should_load_checkpoint)
  {
    load_from_checkpoint = should_load_checkpoint;
  }

  /**
   * @brief Return if checkpoint should be loaded.
   */
  [[nodiscard]] bool
  should_load_checkpoint() const
  {
    return load_from_checkpoint;
  }

  /**
   * @brief Return if the increment should be checkpointted.
   */
  [[nodiscard]] bool
  should_checkpoint(unsigned int increment) const;

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
   * @brief Get the file name.
   */
  [[nodiscard]] const std::string &
  get_file_name() const
  {
    return file_name;
  }

  /**
   * @brief Set the file name
   */
  void
  set_file_name(const std::string &_file_name)
  {
    file_name = _file_name;
  }

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
   * @brief Whether to print timing information with checkpoint
   */
  void
  set_print_timing_with_checkpoint(const bool &_print_timing_with_checkpoint)
  {
    print_timing_with_checkpoint = _print_timing_with_checkpoint;
  }

private:
  // Whether to load from a checkpoint
  bool load_from_checkpoint = false;

  // Checkpoint file name
  std::string file_name;

  // Whether to print timing information with checkpoint
  // TODO (landinjm): Implement this.
  bool print_timing_with_checkpoint = false;

  // List of increments that checkpoint the solution to file
  std::set<unsigned int> checkpoint_list;
};

inline bool
CheckpointParameters::should_checkpoint(unsigned int increment) const
{
  return checkpoint_list.contains(increment);
}

inline void
CheckpointParameters::postprocess_and_validate()
{}

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

PRISMS_PF_END_NAMESPACE
