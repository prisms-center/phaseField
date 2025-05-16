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
 * \brief Struct that holds output parameters.
 */
struct outputParameters
{
public:
  /**
   * \brief Return is the current increment should be output.
   */
  [[nodiscard]] bool
  should_output(unsigned int increment) const;

  /**
   * \brief Postprocess and validate parameters.
   */
  void
  postprocess_and_validate(const temporalDiscretization &temporal_discretization);

  /**
   * \brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  /**
   * \brief Get the file type.
   */
  [[nodiscard]] const std::string &
  get_file_type() const
  {
    return file_type;
  }

  /**
   * \brief Set the file type
   */
  void
  set_file_type(const std::string &_file_type)
  {
    file_type = _file_type;
  }

  /**
   * \brief Get the file name.
   */
  [[nodiscard]] const std::string &
  get_file_name() const
  {
    return file_name;
  }

  /**
   * \brief Set the file name
   */
  void
  set_file_name(const std::string &_file_name)
  {
    file_name = _file_name;
  }

  /**
   * \brief Set the print output period
   */
  void
  set_print_output_period(const unsigned int &_print_output_period)
  {
    print_output_period = _print_output_period;
  }

  /**
   * \brief Get the patch subdivisions
   */
  [[nodiscard]] unsigned int
  get_patch_subdivisions() const
  {
    return patch_subdivisions;
  }

  /**
   * \brief Set the patch subdivisions
   */
  void
  set_patch_subdivisions(const unsigned int &_patch_subdivisions)
  {
    patch_subdivisions = _patch_subdivisions;
  }

  /**
   * \brief Set the output condition
   */
  void
  set_output_condition(const std::string &_condition)
  {
    condition = _condition;
  }

  /**
   * \brief Set the number of outputs
   */
  void
  set_n_outputs(const unsigned int &_n_outputs)
  {
    n_outputs = _n_outputs;
  }

  /**
   * \brief Set the user output list
   */
  void
  set_user_output_list(const std::vector<int> &_user_output_list)
  {
    user_output_list = _user_output_list;
  }

  /**
   * \brief Whether to print timing information with output
   */
  void
  set_print_timing_with_output(const bool &_print_timing_with_output)
  {
    print_timing_with_output = _print_timing_with_output;
  }

private:
  // Output file type
  std::string file_type;

  // Output file name
  std::string file_name;

  // The number of steps between outputting relevant information to screen
  unsigned int print_output_period = UINT_MAX;

  // The number of subdivisions to apply when building patches. By default this is the
  // element degree.
  unsigned int patch_subdivisions = 0;

  // Output condition type
  std::string condition;

  // Number of outputs
  unsigned int n_outputs = 0;

  // User given output_list
  std::vector<int> user_output_list;

  // Whether to print timing information with output
  // TODO (landinjm): Implement this.
  bool print_timing_with_output = false;

  // List of increments that output the solution to file
  std::set<unsigned int> output_list;
};

inline bool
outputParameters::should_output(unsigned int increment) const
{
  return output_list.contains(increment);
}

inline void
outputParameters::postprocess_and_validate(
  const temporalDiscretization &temporal_discretization)
{
  // If the user has specified a list and we have list output use that and return early
  if (condition == "LIST")
    {
      for (const auto &increment : user_output_list)
        {
          output_list.insert(static_cast<unsigned int>(increment));
        }
      return;
    }

  // If the number of outputs is 0 return early
  if (n_outputs == 0)
    {
      return;
    }

  // If the number of outputs is greater than the number of increments, force them to be
  // equivalent
  n_outputs = std::min(n_outputs, temporal_discretization.get_total_increments());

  // Determine the output list from the other criteria
  if (condition == "EQUAL_SPACING")
    {
      for (unsigned int iteration = 0;
           iteration <= temporal_discretization.get_total_increments();
           iteration += temporal_discretization.get_total_increments() / n_outputs)
        {
          output_list.insert(iteration);
        }
    }
  else if (condition == "LOG_SPACING")
    {
      output_list.insert(0);
      for (unsigned int output = 1; output <= n_outputs; output++)
        {
          output_list.insert(static_cast<unsigned int>(std::round(
            std::pow(static_cast<double>(temporal_discretization.get_total_increments()),
                     static_cast<double>(output) / static_cast<double>(n_outputs)))));
        }
    }
  else if (condition == "N_PER_DECADE")
    {
      AssertThrow(temporal_discretization.get_total_increments() > 1,
                  dealii::ExcMessage("For n per decaded spaced outputs, the number of "
                                     "increments must be greater than 1."));

      output_list.insert(0);
      output_list.insert(1);
      for (unsigned int iteration = 2;
           iteration <= temporal_discretization.get_total_increments();
           iteration++)
        {
          const auto decade = static_cast<unsigned int>(std::ceil(std::log10(iteration)));
          const auto step_size =
            static_cast<unsigned int>(std::pow(10, decade) / n_outputs);
          if (iteration % step_size == 0)
            {
              output_list.insert(iteration);
            }
        }
    }
  else
    {
      AssertThrow(false, UnreachableCode());
    }
}

inline void
outputParameters::print_parameter_summary() const
{
  conditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Output Parameters\n"
    << "================================================\n"
    << "Output file type: " << file_type << "\n"
    << "Output file name: " << file_name << "\n"
    << "Output subdivisions: " << patch_subdivisions << "\n"
    << "Print output period: " << print_output_period << "\n"
    << "Output condition: " << condition << "\n"
    << "Number of outputs: " << n_outputs << "\n"
    << "Print timing info: " << bool_to_string(print_timing_with_output) << "\n";

  conditionalOStreams::pout_summary() << "Output iteration list: ";
  for (const auto &iteration : output_list)
    {
      conditionalOStreams::pout_summary() << iteration << " ";
    }
  conditionalOStreams::pout_summary() << "\n\n" << std::flush;
}

PRISMS_PF_END_NAMESPACE
