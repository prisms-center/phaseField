// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/utilities.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/parameter_base.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include <climits>
#include <execution>
#include <set>
#include <string>
#include <unordered_map>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Simple struct for field output.
 */
struct FieldOutputParameters
{
  /**
   * @brief VTK output types.
   */
  enum OutputType : std::uint8_t
  {
    VTU,
    VTK,
    PVTU
  };

  /**
   * @brief File type for field output.
   *
   * This determines what type of files are output. For examples, vtu, pvtu, and vtk.
   * The pvtu files are written in parallel and should be the fastest to output when
   * working with large simulations; however, they can become numerous and harder to
   * visualize when you have many threads, as each thread outputs.
   */
  OutputType file_type = OutputType::VTU;

  /**
   * @brief Number of subdivisions to apply when building patches.
   *
   * By default this is the element degree. You can choose higher or lower in order to
   * control how much interpolation should be done to the fields when outputting them.
   * The element degree is the default value to show the increase in DoFs and accuracy
   * when using higher order elements.
   */
  unsigned int patch_subdivisions = 0;

  /**
   * @brief Compression level for output
   *
   * Since zlib is always available, users can select what sort of compression to use for
   * binary output. The options are no compression, best speed, best size, and default (a
   * balance of speed and size).
   */
  dealii::DataOutBase::CompressionLevel compression_level =
    dealii::DataOutBase::CompressionLevel::default_compression;

  /**
   * @brief Folder for field output.
   *
   * Choosing this will give you your results like outputs/solution-*.vtu
   */
  std::string folder;

  /**
   * @brief Base filename for field output.
   *
   * This is the base filename for the outputs. For example, solution-*.vtu or
   * file_name-*.vtu
   */
  std::string file_name;

  /**
   * @brief A list of output steps.
   *
   * This is determined by a combination of the number of outputs and the total number of
   * steps. When we reach a step contained in the list, we output.
   */
  std::set<unsigned int> output_list;

  /**
   * @brief A list of fields that are output.
   *
   * This determines which fields are output during output steps. Each entry corresponds
   * to a field index and the dependency type. Typically, the dependency type will always
   * be normal; however, one may want to visualize change and old fields for debugging.
   *
   * Additionally, this variable is useful when runnings simulation with multiple order
   * parameters. Most times, we don't need to output all of the order parameters because
   * we only care about combined fields (e.g., grain ids). When we don't output all of the
   * order parameters we save lots of disk space because we go from n fields to 1.
   */
  std::set<std::pair<Types::Index, DependencyType>> output_fields;
};

/**
 * @brief A class that determines how often and what fields are output.
 */
class FieldOutputParameterLoader : ParameterBase
{
public:
  /**
   * @brief Declare parameters.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler) const override
  {
    parameter_handler.enter_subsection("field output");

    parameter_handler.declare_entry("folder",
                                    "",
                                    dealii::Patterns::Anything(),
                                    "The folder where output files are stored.");
    parameter_handler.declare_entry("file name",
                                    "solution",
                                    dealii::Patterns::Anything(),
                                    "The base name for the output file, before the time "
                                    "step and processor info are added.");
    parameter_handler.declare_entry("file type",
                                    "vtu",
                                    dealii::Patterns::Selection("vtu|vtk|pvtu"),
                                    "The output file type (either vtu, pvtu, or vtk).");
    parameter_handler.declare_entry(
      "subdivisions",
      "0",
      dealii::Patterns::Integer(0, INT_MAX),
      "The number of subdivisions to apply to the mesh when building output patches. If "
      "0, the degree is used.");
    parameter_handler.declare_entry("compression level",
                                    "default",
                                    dealii::Patterns::Selection(
                                      "default|best speed|best size|none"),
                                    "The compression level of the output (either "
                                    "default, best speed, best size, or none).");

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
      "The list of time steps to output. Used for the LIST type only and must be comma "
      "delimited.");
    parameter_handler.declare_entry("number",
                                    "10",
                                    dealii::Patterns::Integer(0, INT_MAX),
                                    "The number of outputs (or number of outputs "
                                    "per decade for the N_PER_DECADE type).");

    parameter_handler.declare_entry(
      "fields",
      "",
      dealii::Patterns::List(dealii::Patterns::Anything(), 0, INT_MAX, ","),
      "The list of fields to output. Must be comma delimited. Additionally, for the "
      "output of change and old fields, they must follow the same delimiters in "
      "VariableAttributes. In other words, something like `set fields = n1, old_1(n1), "
      "change(n1)`.");

    parameter_handler.leave_subsection();
  }

  /**
   * @brief Read parameters.
   */
  void
  read(dealii::ParameterHandler &parameter_handler) override
  {
    const std::unordered_map<std::string, FieldOutputParameters::OutputType>
      output_type_table = {
        {"vtu",  FieldOutputParameters::OutputType::VTU },
        {"vtk",  FieldOutputParameters::OutputType::VTK },
        {"pvtu", FieldOutputParameters::OutputType::PVTU}
    };
    const std::unordered_map<std::string, dealii::DataOutBase::CompressionLevel>
      compression_level_table = {
        {"default",    dealii::DataOutBase::CompressionLevel::default_compression},
        {"best speed", dealii::DataOutBase::CompressionLevel::best_speed         },
        {"best size",  dealii::DataOutBase::CompressionLevel::best_compression   },
        {"none",       dealii::DataOutBase::CompressionLevel::no_compression     }
    };

    parameter_handler.enter_subsection("field output");

    parameters.folder    = parameter_handler.get("folder");
    parameters.file_name = parameter_handler.get("file name");
    parameters.file_type =
      string_to_type(parameter_handler.get("file type"), output_type_table);
    parameters.patch_subdivisions = parameter_handler.get_integer("subdivisions");
    parameters.compression_level =
      string_to_type(parameter_handler.get("compression level"), compression_level_table);

    condition        = parameter_handler.get("condition");
    user_output_list = dealii::Utilities::string_to_int(
      dealii::Utilities::split_string_list(parameter_handler.get("list")));
    n_outputs  = parameter_handler.get_integer("number");
    field_list = dealii::Utilities::split_string_list(parameter_handler.get("fields"));

    parameter_handler.leave_subsection();
  }

  /**
   * @brief Postprocess parameters.
   */
  void
  postprocess(unsigned int n_increments)
  {
    fill_out_output_list(n_increments);

    // Determine which fields we are outputting.
  }

  /**
   * @brief Validate parameters.
   */
  void
  validate(const std::vector<VariableAttributes> &field_attributes) const override
  {}

  /**
   * @brief Package parameters.
   */
  [[nodiscard]] FieldOutputParameters
  package() const
  {
    return parameters;
  }

  /**
   * @brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const override
  {
    ConditionalOStreams::pout_summary()
      << "================================================\n"
      << "  Field Output Parameter Class\n"
      << "================================================\n"
      << std::flush;
  }

private:
  void
  fill_out_output_list(unsigned int n_increments)
  {
    // If the user has specified a list and we have list output use that and return early
    if (condition == "LIST")
      {
        for (const auto &increment : user_output_list)
          {
            parameters.output_list.insert(static_cast<unsigned int>(increment));
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
    n_outputs = std::min(n_outputs, n_increments);

    // If the number of increments is 0, we only output the initial condition. We can set
    // that and return early.
    if (n_increments == 0)
      {
        parameters.output_list.insert(0);
        return;
      }

    // Determine the output list from the other criteria
    if (condition == "EQUAL_SPACING")
      {
        for (unsigned int iteration = 0; iteration <= n_increments;
             iteration += n_increments / n_outputs)
          {
            parameters.output_list.insert(iteration);
          }
      }
    else if (condition == "LOG_SPACING")
      {
        parameters.output_list.insert(0);
        for (unsigned int output = 1; output <= n_outputs; output++)
          {
            parameters.output_list.insert(static_cast<unsigned int>(std::round(
              std::pow(static_cast<double>(n_increments),
                       static_cast<double>(output) / static_cast<double>(n_outputs)))));
          }
      }
    else if (condition == "N_PER_DECADE")
      {
        AssertThrow(n_increments > 1,
                    dealii::ExcMessage("For n per decaded spaced outputs, the number of "
                                       "increments must be greater than 1."));
        parameters.output_list.insert(0);
        parameters.output_list.insert(1);
        for (unsigned int iteration = 2; iteration <= n_increments; iteration++)
          {
            const auto decade =
              static_cast<unsigned int>(std::ceil(std::log10(iteration)));
            const auto step_size =
              static_cast<unsigned int>(std::pow(10, decade) / n_outputs);
            if (iteration % step_size == 0)
              {
                parameters.output_list.insert(iteration);
              }
          }
      }
    else
      {
        AssertThrow(false, UnreachableCode());
      }
  }

  FieldOutputParameters parameters;

  std::string condition;

  std::vector<int> user_output_list;

  unsigned int n_outputs = 0;

  std::vector<std::string> field_list;
};

PRISMS_PF_END_NAMESPACE
