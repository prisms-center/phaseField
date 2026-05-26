// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
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
 * @brief Struct that holds output parameters.
 */
struct OutputParameters
{
  /**
   * @brief Return if the increment should be outputted.
   */
  [[nodiscard]] bool
  should_output(unsigned int increment) const;

  /**
   * @brief Postprocess and validate parameters.
   */
  void
  validate();

  /**
   * @brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

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
    if (num_increments == 0)
      {
        output_list.insert(0);
        return;
      }
    if (!num_outputs)
      {
        return;
      }
    unsigned int period = std::max(1U, num_increments / num_outputs);
    for (unsigned int output = 0; output <= num_increments; output += period)
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
    if (!num_outputs)
      {
        return;
      }
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
    if (!num_outputs)
      {
        return;
      }
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

  // Output directory
  std::string directory = "solutions";

  // Output file type ()
  std::string file_type = "vtu";

  // Output file name
  std::string file_name = "solution";

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
OutputParameters::validate()
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

inline void
OutputParameters::declare_parameters(dealii::ParameterHandler &parameter_handler) const
{
  parameter_handler.enter_subsection("output");
  {
    parameter_handler.declare_entry("directory",
                                    "solutions",
                                    dealii::Patterns::Anything(),
                                    "The name of the output directory.");
    parameter_handler.declare_entry("file name",
                                    "solution",
                                    dealii::Patterns::Anything(),
                                    "The prefix of the output files, before the "
                                    "time step and processor info are added.");
    parameter_handler.declare_entry(
      "file type",
      "vtu",
      dealii::Patterns::Selection("vtu|vtk|pvtu|xdmf"),
      "The output file type (either vtu, pvtu, vtk, or xdmf).");
    parameter_handler.declare_entry(
      "subdivisions",
      "0",
      dealii::Patterns::Integer(0, INT_MAX),
      "The number of subdivisions to apply to the mesh when building output patches.");
    parameter_handler.declare_entry(
      "condition",
      "EQUAL_SPACING",
      dealii::Patterns::Selection("EQUAL_SPACING|LOG_SPACING|N_PER_DECADE|LIST"),
      "The spacing type for outputting the solution fields (either EQUAL_SPACING, "
      "LOG_SPACING, N_PER_DECADE, or LIST).");
    parameter_handler.declare_entry(
      "list",
      "0",
      dealii::Patterns::List(dealii::Patterns::Integer(0, INT_MAX), 0, INT_MAX, ","),
      "The list of time steps to output, used for the LIST type.");
    parameter_handler.declare_entry("number",
                                    "10",
                                    dealii::Patterns::Integer(0, INT_MAX),
                                    "The number of outputs (or number of outputs "
                                    "per decade for the N_PER_DECADE type).");
    parameter_handler.declare_entry(
      "print step period",
      "2147483647",
      dealii::Patterns::Integer(1, INT_MAX),
      "The number of time steps between updates to the screen.");
    parameter_handler.declare_entry(
      "timing information with output",
      "false",
      dealii::Patterns::Bool(),
      "Whether to print the summary table of the wall time and wall time for "
      "individual subroutines every time the code outputs.");
    parameter_handler.declare_alias("directory", "folder name");
  }
  parameter_handler.leave_subsection();
}

inline void
OutputParameters::assign_parameters(dealii::ParameterHandler     &parameter_handler,
                                    const TemporalDiscretization &temporal_discretization)
{
  parameter_handler.enter_subsection("output");
  {
    directory = parameter_handler.get("directory");
    file_name = parameter_handler.get("file name");
    file_type = parameter_handler.get("file type");
    patch_subdivisions =
      static_cast<unsigned int>(parameter_handler.get_integer("subdivisions"));

    std::string  condition = parameter_handler.get("condition");
    unsigned int num_outputs =
      static_cast<unsigned int>(parameter_handler.get_integer("number"));
    unsigned int num_increments = temporal_discretization.num_increments;
    if (condition == "EQUAL_SPACING")
      {
        add_equal_spacing_outputs(num_outputs, num_increments);
      }
    else if (condition == "LOG_SPACING")
      {
        add_log_spacing_outputs(num_outputs, num_increments);
      }
    else if (condition == "N_PER_DECADE")
      {
        add_n_per_decade_outputs(num_outputs, num_increments);
      }
    add_output_list(dealii::Utilities::string_to_int(
      dealii::Utilities::split_string_list(parameter_handler.get("list"))));

    print_output_period =
      static_cast<unsigned int>(parameter_handler.get_integer("print step period"));
    print_timing_with_output =
      parameter_handler.get_bool("timing information with output");
  }
  parameter_handler.leave_subsection();
}

PRISMS_PF_END_NAMESPACE
