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
 * @brief Struct that holds output parameters.
 */
struct OutputParameters
{
public:
  /**
   * @brief Return if the increment should be outputted.
   */
  [[nodiscard]] bool
  should_output(unsigned int increment) const;

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
   * @brief Get the file type.
   */
  [[nodiscard]] const std::string &
  get_file_type() const
  {
    return file_type;
  }

  /**
   * @brief Set the file type
   */
  void
  set_file_type(const std::string &_file_type)
  {
    file_type = _file_type;
  }

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
   * @brief Set the print output period
   */
  void
  set_print_output_period(const unsigned int &_print_output_period)
  {
    print_output_period = _print_output_period;
  }

  /**
   * @brief Get the patch subdivisions
   */
  [[nodiscard]] unsigned int
  get_patch_subdivisions() const
  {
    return patch_subdivisions;
  }

  /**
   * @brief Set the patch subdivisions
   */
  void
  set_patch_subdivisions(const unsigned int &_patch_subdivisions)
  {
    patch_subdivisions = _patch_subdivisions;
  }

  /**
   * @brief Set the user output list
   */
  template <typename ListType>
  void
  add_output_list(const ListType &list)
  {
    output_list.insert(list.begin(), list.end());
  }

  /**
   * @brief Set the user output list
   */
  void
  add_equal_spacing_outputs(unsigned int num_outputs, unsigned int num_increments)
  {
    for (unsigned int output = 0; output <= num_increments;
         output += num_increments / num_outputs)
      {
        output_list.insert(output);
      }
  }

  /**
   * @brief Set the user output list
   */
  void
  add_log_spacing_outputs(unsigned int num_outputs, unsigned int num_increments)
  {
    for (unsigned int output = 1; output <= num_outputs; output++)
      {
        output_list.insert(static_cast<unsigned int>(std::round(
          std::pow(double(num_increments), double(output) / double(num_outputs)))));
      }
  }

  /**
   * @brief Set the user output list
   */
  void
  add_n_per_decade_outputs(unsigned int num_outputs, unsigned int num_increments)
  {
    AssertThrow(num_increments > 1,
                dealii::ExcMessage("For n per decaded spaced outputs, the number of "
                                   "increments must be greater than 1."));

    output_list.insert(0);
    output_list.insert(1);
    for (unsigned int iteration = 2; iteration <= num_increments; iteration++)
      {
        const auto decade = static_cast<unsigned int>(std::ceil(std::log10(iteration)));
        const auto step_size =
          static_cast<unsigned int>(std::pow(10, decade) / num_outputs);
        if (iteration % step_size == 0)
          {
            output_list.insert(iteration);
          }
      }
  }

  /**
   * @brief Set the user output list
   */
  void
  clear_output_list()
  {
    output_list = {0};
  }

  /**
   * @brief Get the number of outputs that will be made
   */
  [[nodiscard]] unsigned int
  get_num_outputs() const
  {
    return output_list.size();
  }

  /**
   * @brief Whether to print timing information with output
   */
  void
  set_print_timing_with_output(const bool &_print_timing_with_output)
  {
    print_timing_with_output = _print_timing_with_output;
  }

private:
  // Output file type ()
  std::string file_type;

  // Output file name
  std::string file_name;

  // The number of subdivisions to apply when building patches. By default this is the
  // element degree.
  unsigned int patch_subdivisions = 0;

  // The number of steps between outputting relevant information to screen
  unsigned int print_output_period = UINT_MAX;

  // Whether to print timing information with output
  // TODO (landinjm): Implement this.
  bool print_timing_with_output = false;

  // List of increments that output the solution to file
  std::set<unsigned int> output_list = {0};
};

inline bool
OutputParameters::should_output(unsigned int increment) const
{
  return output_list.contains(increment);
}

inline void
OutputParameters::postprocess_and_validate()
{}

inline void
OutputParameters::print_parameter_summary() const
{
  ConditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Output Parameters\n"
    << "================================================\n"
    << "Output file type: " << file_type << "\n"
    << "Output file name: " << file_name << "\n"
    << "Output subdivisions: " << patch_subdivisions << "\n"
    << "Print output period: " << print_output_period << "\n"
    << "Number of outputs: " << get_num_outputs() << "\n"
    << "Print timing info: " << bool_to_string(print_timing_with_output) << "\n";

  ConditionalOStreams::pout_summary() << "Output increment list: ";
  for (const auto &iteration : output_list)
    {
      ConditionalOStreams::pout_summary() << iteration << " ";
    }
  ConditionalOStreams::pout_summary() << "\n\n" << std::flush;
}

PRISMS_PF_END_NAMESPACE
