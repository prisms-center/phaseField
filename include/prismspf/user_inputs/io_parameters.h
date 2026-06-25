// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/utilities.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include <climits>
#include <concepts>
#include <execution>
#include <set>
#include <string>
#include <unordered_map>

PRISMS_PF_BEGIN_NAMESPACE

inline void
add_equal_spacing_outputs(unsigned int            n_outputs,
                          unsigned int            n_increments,
                          std::set<unsigned int> &output_list)
{
  output_list.clear();
  if (n_increments == 0)
    {
      output_list.insert(0);
      return;
    }
  if (!n_outputs)
    {
      return;
    }

  unsigned int period = std::max(1U, n_increments / n_outputs);
  for (unsigned int output = 0; output <= n_increments; output += period)
    {
      output_list.insert(output);
    }
}

inline void
add_log_spacing_outputs(unsigned int            n_outputs,
                        unsigned int            n_increments,
                        std::set<unsigned int> &output_list)
{
  output_list.clear();
  if (n_increments == 0)
    {
      output_list.insert(0);
      return;
    }
  if (!n_outputs)
    {
      return;
    }

  for (unsigned int output = 1; output <= n_outputs; output++)
    {
      output_list.insert((unsigned int) (std::round(
        std::pow(double(n_increments), double(output) / double(n_outputs)))));
    }
}

inline void
add_n_per_decade_outputs(unsigned int            n_outputs,
                         unsigned int            n_increments,
                         std::set<unsigned int> &output_list)
{
  output_list.clear();
  if (n_increments == 0)
    {
      output_list.insert(0);
      return;
    }
  if (!n_outputs)
    {
      return;
    }

  output_list.insert(0);
  output_list.insert(1);
  for (unsigned int iteration = 2; iteration <= n_increments; iteration++)
    {
      const auto decade    = (unsigned int) (std::ceil(std::log10(iteration)));
      const auto step_size = (unsigned int) (std::pow(10, decade) / n_outputs);
      if (iteration % step_size == 0)
        {
          output_list.insert(iteration);
        }
    }
}

template <typename T1, typename T2>
requires std::convertible_to<T1, T2>
inline void
add_list_outputs(const std::vector<T1> &list, std::set<T2> &output_list)
{
  for (const auto &val : list)
    {
      output_list.insert(T2(val));
    }
}

/**
 * @brief Simple struct for field output.
 */
struct FieldOutputParameters : public ParameterBase
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int max_criteria = Numbers::max_subsections) const override
  {
    parameter_handler.enter_subsection("output");
    {
      parameter_handler.declare_entry(
        "file type",
        "vtu",
        dealii::Patterns::Selection("vtu|vtk|pvtu|xdmf"),
        "The output file type (either vtu, pvtu, vtk, or xdmf).");

      parameter_handler.declare_entry(
        "subdivisions",
        "0",
        dealii::Patterns::Integer(0, INT_MAX),
        "The number of subdivisions to apply to the mesh when building output patches. "
        "If 0, the degree is used.");

      parameter_handler.declare_entry("compression level",
                                      "default",
                                      dealii::Patterns::Selection(
                                        "default|speed|size|none"),
                                      "The compression level for output.");

      parameter_handler.declare_entry("directory",
                                      "outputs",
                                      dealii::Patterns::Anything(),
                                      "The name of the output directory.");
      parameter_handler.declare_alias("directory", "folder name");

      parameter_handler.declare_entry("file name",
                                      "solution",
                                      dealii::Patterns::Anything(),
                                      "The prefix of the output files, before the "
                                      "time step and processor info are added.");

      parameter_handler.declare_entry(
        "condition",
        "EQUAL_SPACING",
        dealii::Patterns::Selection("EQUAL_SPACING|LOG_SPACING|N_PER_DECADE|LIST"),
        "The spacing type for outputting the solution fields.");
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
        "variables",
        "",
        dealii::Patterns::List(dealii::Patterns::Anything(), 0, INT_MAX, ","),
        "The list of the fields to output. Must be comma delimited. Additionally, for "
        "the output of left-hand side and old fields, they must follow the same "
        "delimiters that are used in dependency sets. In other words, something like "
        "`set variables = n1, old_1(n1), lhs(n1)`.");
    }
    parameter_handler.leave_subsection();
  };

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              n_increments,
         unsigned int              max_criteria = Numbers::max_subsections)
  {
    const static std::unordered_map<std::string, FieldOutputParameters::OutputType>
      output_type_table = {
        {"vtu",  FieldOutputParameters::OutputType::VTU },
        {"vtk",  FieldOutputParameters::OutputType::VTK },
        {"pvtu", FieldOutputParameters::OutputType::PVTU},
        {"xdmf", FieldOutputParameters::OutputType::XDMF}
    };
    const static std::unordered_map<std::string, dealii::DataOutBase::CompressionLevel>
      compression_level_table = {
        {"default",    dealii::DataOutBase::CompressionLevel::default_compression},
        {"best speed", dealii::DataOutBase::CompressionLevel::best_speed         },
        {"best size",  dealii::DataOutBase::CompressionLevel::best_compression   },
        {"none",       dealii::DataOutBase::CompressionLevel::no_compression     }
    };

    parameter_handler.enter_subsection("output");
    {
      folder             = parameter_handler.get("directory");
      file_name          = parameter_handler.get("file name");
      patch_subdivisions = (unsigned int) (parameter_handler.get_integer("subdivisions"));

      file_type = output_type_table.at(parameter_handler.get("file type"));

      compression_level =
        compression_level_table.at(parameter_handler.get("compression level"));

      add_list_outputs(dealii::Utilities::split_string_list(
                         parameter_handler.get("variables")),
                       output_fields);

      std::string  condition = parameter_handler.get("condition");
      unsigned int n_outputs = (unsigned int) (parameter_handler.get_integer("number"));

      if (condition == "EQUAL_SPACING")
        {
          add_equal_spacing_outputs(n_outputs, n_increments, output_list);
        }
      else if (condition == "LOG_SPACING")
        {
          add_log_spacing_outputs(n_outputs, n_increments, output_list);
        }
      else if (condition == "N_PER_DECADE")
        {
          add_n_per_decade_outputs(n_outputs, n_increments, output_list);
        }
      else
        {
          add_list_outputs(dealii::Utilities::string_to_int(
                             dealii::Utilities::split_string_list(
                               parameter_handler.get("list"))),
                           output_list);
        }
    }
    parameter_handler.leave_subsection();
  };

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override {
    // TODO: Do this later
  };

  /**
   * @brief VTK output types.
   */
  enum OutputType : std::uint8_t
  {
    VTU,
    VTK,
    PVTU,
    XDMF
  };

  /**
   * @brief Whether a given increment should be outputted.
   */
  [[nodiscard]] bool
  should_output(unsigned int increment) const
  {
    return output_list.contains(increment);
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
  std::string folder = "outputs";

  /**
   * @brief Base filename for field output.
   *
   * This is the base filename for the outputs. For example, solution-*.vtu or
   * file_name-*.vtu
   */
  std::string file_name = "solution";

  /**
   * @brief A list of output steps.
   *
   * This is determined by a combination of the number of outputs and the total number of
   * steps. When we reach a step contained in the list, we output.
   */
  std::set<unsigned int> output_list = {0};

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
  std::set<std::string> output_fields;
};

PRISMS_PF_END_NAMESPACE
