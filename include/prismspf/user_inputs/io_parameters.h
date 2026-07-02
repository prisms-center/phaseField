// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/utilities.h>

#include <boost/serialization/vector.hpp>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/parameter_base.h>

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
          unsigned int max_criteria = Numbers::max_subsections) const override;

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              n_increments,
         unsigned int              max_criteria = Numbers::max_subsections);

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override;

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
  should_output(unsigned int increment) const;

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

/**
 * @brief Simple struct for restart output.
 */
struct RestartOutputParameters : public ParameterBase
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int max_criteria = Numbers::max_subsections) const override;

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              n_increments,
         unsigned int              max_criteria = Numbers::max_subsections);

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override;

  /**
   * @brief Whether a given increment should be outputted.
   */
  [[nodiscard]] bool
  should_output(unsigned int increment) const;

  /**
   * @brief Whether to load from a checkpoint
   *
   * TODO: This doesn't really belong here. Maybe it does? With loading from restart it
   * overrides any other input conditions so it might not fit super well in the input
   * section either. I don't see a good place to put it and it somewhat fits here so
   * leaving it.
   */
  bool load_from_checkpoint = false;

  /**
   * @brief Folder for checkpoint output.
   *
   * Choosing this will give you your results like outputs/solution-*.vtu
   */
  std::string folder = "outputs";

  /**
   * @brief Base filename for checkpoint output.
   *
   * This is the base filename for the outputs. For example, solution-*.vtu or
   * file_name-*.vtu
   */
  std::string file_name = "checkpoint";

  /**
   * @brief A list of output steps.
   *
   * This is determined by a combination of the number of outputs and the total number of
   * steps. When we reach a step contained in the list, we output.
   */
  std::set<unsigned int> output_list = {0};
};

/**
 * @brief Initial condition file
 */
struct InitialConditionFile
{
  /**
   * @brief Data formats for input initial conditions.
   */
  enum DataFormatType : std::uint8_t
  {
    FlatBinary,
    VTKUnstructuredGrid,
    VTKXMLUnstructuredGrid,
    VTKPXMLUnstructuredGrid,
    VTKXMLImageData
  };

  // File name
  std::string file_name;

  // Data format
  DataFormatType format;

  // File variable names
  std::vector<std::string> file_variable_names;

  // Simulation variable names
  std::vector<std::string> simulation_variable_names;

  // Number of data points in each direction
  std::array<unsigned int, 3> n_data_points = {
    {0, 0, 0}
  };
};

/**
 * @brief Simple struct for field input.
 */
struct FieldInputParameters : public ParameterBase
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int max_criteria = Numbers::max_subsections) const override;

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              max_criteria = Numbers::max_subsections) override;

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override;

  /**
   * @brief Whether to load initial conditions from file
   */
  bool load_from_file = false;

  // Collection of initial condition files
  std::vector<InitialConditionFile> initial_condition_files;
};

PRISMS_PF_END_NAMESPACE
