// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/data_out.h>

#include <prismspf/core/solution_output.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <fstream>
#include <iomanip>
#include <map>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, typename number>
SolutionOutput<dim, number>::SolutionOutput(const VectorType               &solution,
                                            const dealii::DoFHandler<dim>  &dof_handler,
                                            const unsigned int             &degree,
                                            const std::string              &name,
                                            const UserInputParameters<dim> &user_inputs)
{
  // Some stuff to determine the actual name of the output file.
  const auto n_trailing_digits = static_cast<unsigned int>(
    std::floor(
      std::log10(user_inputs.get_temporal_discretization().get_total_increments())) +
    1);

  // Init data out
  dealii::DataOut<dim> data_out;

  solution.update_ghost_values();

  // Add data vector
  data_out.add_data_vector(dof_handler, solution, name);

  solution.zero_out_ghost_values();

  // Build patches to linearly interpolate from higher order element degrees. Note that
  // this essentially converts the element to an equal amount of subdivisions in the
  // output. This does not make subdivisions and element degree equivalent in the
  // simulation!
  const unsigned int n_divisions =
    user_inputs.get_output_parameters().get_patch_subdivisions() == 0
      ? degree
      : user_inputs.get_output_parameters().get_patch_subdivisions();
  data_out.build_patches(n_divisions);

  // Set some flags for data output
  dealii::DataOutBase::VtkFlags flags;
  flags.time                = user_inputs.get_temporal_discretization().get_time();
  flags.cycle               = user_inputs.get_temporal_discretization().get_increment();
  flags.print_date_and_time = true;
#ifdef PRISMS_PF_WITH_ZLIB
  // TODO (landinjm): Make this a user input parameter so they can select between
  // compression levels. Best speed can be the default
  flags.compression_level = dealii::DataOutBase::CompressionLevel::best_speed;
#endif
  data_out.set_flags(flags);

  // Write to file based on the user input.
  const std::string  directory = "./";
  const unsigned int increment =
    user_inputs.get_temporal_discretization().get_increment();

  if (user_inputs.get_output_parameters().get_file_type() == "vtu")
    {
      std::ostringstream increment_stream;
      increment_stream << std::setw(static_cast<int>(n_trailing_digits))
                       << std::setfill('0') << increment;
      const std::string filename =
        directory + name + "_" + increment_stream.str() + ".vtu";
      data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD);
    }
  else if (user_inputs.get_output_parameters().get_file_type() == "pvtu")
    {
      data_out.write_vtu_with_pvtu_record(directory,
                                          name,
                                          increment,
                                          MPI_COMM_WORLD,
                                          n_trailing_digits);
    }
  else if (user_inputs.get_output_parameters().get_file_type() == "vtk")
    {
      std::ostringstream increment_stream;
      increment_stream << std::setw(static_cast<int>(n_trailing_digits))
                       << std::setfill('0') << increment;
      const std::string filename =
        directory + name + "_" + increment_stream.str() + ".vtk";
      std::ofstream vtk_output(filename);
      data_out.write_vtk(vtk_output);
    }
  else
    {
      AssertThrow(false, UnreachableCode());
    }

  // Update the ghost values again to allow for read access
  solution.update_ghost_values();
}

template <unsigned int dim, typename number>
SolutionOutput<dim, number>::SolutionOutput(
  const std::map<unsigned int, VectorType *>         &solution_set,
  const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers,
  const unsigned int                                 &degree,
  const std::string                                  &name,
  const UserInputParameters<dim>                     &user_inputs)
{
  // Some stuff to determine the actual name of the output file.
  const auto n_trailing_digits = static_cast<unsigned int>(
    std::floor(
      std::log10(user_inputs.get_temporal_discretization().get_total_increments())) +
    1);

  // Init data out
  dealii::DataOut<dim> data_out;

  // Add data vectors
  for (const auto &[index, variable] : user_inputs.get_variable_attributes())
    {
      auto *solution = solution_set.at(index);
      solution->update_ghost_values();

      // Mark field as Scalar/Vector
      const bool is_scalar =
        variable.field_info.tensor_rank == FieldInfo::TensorRank::Scalar;
      const unsigned int n_components = is_scalar ? 1 : dim;

      const std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
        data_type(n_components,
                  is_scalar
                    ? dealii::DataComponentInterpretation::component_is_scalar
                    : dealii::DataComponentInterpretation::component_is_part_of_vector);

      const std::vector<std::string> names(n_components, variable.get_name());

      data_out.add_data_vector(*(dof_handlers.at(index)), *solution, names, data_type);

      solution->zero_out_ghost_values();
    }

  // Build patches to linearly interpolate from higher order element degrees. Note that
  // this essentially converts the element to an equal amount of subdivisions in the
  // output. This does not make subdivisions and element degree equivalent in the
  // simulation!
  const unsigned int n_divisions =
    user_inputs.get_output_parameters().get_patch_subdivisions() == 0
      ? degree
      : user_inputs.get_output_parameters().get_patch_subdivisions();
  data_out.build_patches(n_divisions);

  // Set some flags for data output
  dealii::DataOutBase::VtkFlags flags;
  flags.time                = user_inputs.get_temporal_discretization().get_time();
  flags.cycle               = user_inputs.get_temporal_discretization().get_increment();
  flags.print_date_and_time = true;
#ifdef PRISMS_PF_WITH_ZLIB
  flags.compression_level = dealii::DataOutBase::CompressionLevel::best_speed;
#endif
  data_out.set_flags(flags);

  // Write to file based on the user input.
  const std::string  directory = "./";
  const unsigned int increment =
    user_inputs.get_temporal_discretization().get_increment();

  if (user_inputs.get_output_parameters().get_file_type() == "vtu")
    {
      std::ostringstream increment_stream;
      increment_stream << std::setw(static_cast<int>(n_trailing_digits))
                       << std::setfill('0') << increment;
      const std::string filename =
        directory + name + "_" + increment_stream.str() + ".vtu";
      data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD);
    }
  else if (user_inputs.get_output_parameters().get_file_type() == "pvtu")
    {
      data_out.write_vtu_with_pvtu_record(directory,
                                          name,
                                          increment,
                                          MPI_COMM_WORLD,
                                          n_trailing_digits);
    }
  else if (user_inputs.get_output_parameters().get_file_type() == "vtk")
    {
      std::ostringstream increment_stream;
      increment_stream << std::setw(static_cast<int>(n_trailing_digits))
                       << std::setfill('0') << increment;
      const std::string filename =
        directory + name + "_" + increment_stream.str() + ".vtk";
      std::ofstream vtk_output(filename);
      data_out.write_vtk(vtk_output);
    }
  else
    {
      AssertThrow(false, UnreachableCode());
    }

  // Update the ghost values again to allow for read access
  for (const auto &[index, variable] : user_inputs.get_variable_attributes())
    {
      auto *solution = solution_set.at(index);
      solution->update_ghost_values();
    }
}

#include "core/solution_output.inst"

PRISMS_PF_END_NAMESPACE
