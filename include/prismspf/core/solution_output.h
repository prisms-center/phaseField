// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/data_out.h>

#include <prismspf/core/dof_manager.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/simulation_timer.h>
#include <prismspf/core/solution_indexer.h>
#include <prismspf/core/solution_output.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/solvers/solve_context.h>

#include <prismspf/user_inputs/output_parameters.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

/**
 * @brief Class that outputs a passed solution to vtu, vtk, or pvtu
 */
template <unsigned int dim, typename number>
class SolutionOutput
{
public:
  using VectorType = dealii::LinearAlgebra::distributed::Vector<number>;

  /**
   * @brief Outputs all fields in the solution set.
   */
  SolutionOutput(const std::vector<FieldAttributes> &field_attributes,
                 const SolutionIndexer<dim, number> &solution_indexer,
                 const SimulationTimer              &sim_timer,
                 const DofManager<dim>              &dof_manager,
                 const unsigned int                 &degree,
                 const std::string                  &file_prefix,
                 const UserInputParameters<dim>     &user_inputs)
  {
    const OutputParameters &output_parameters = user_inputs.get_output_parameters();
    // Some stuff to determine the actual name of the output file.
    const auto n_trailing_digits = static_cast<unsigned int>(
      std::floor(std::log10(user_inputs.get_temporal_discretization().num_increments)) +
      1);

    // Init data out
    dealii::DataOut<dim> data_out;

    // Add data vectors
    for (unsigned int index = 0; index < field_attributes.size(); ++index)
      {
        const FieldAttributes &field    = field_attributes[index];
        auto                  &solution = solution_indexer.get_solution_vector(index);
        solution.update_ghost_values();

        // Mark field as Scalar/Vector
        const unsigned int n_components =
          (field.field_type == TensorRank::Scalar) ? 1 : dim;

        const std::vector<
          dealii::DataComponentInterpretation::DataComponentInterpretation>
          data_type(n_components,
                    (n_components == 1)
                      ? dealii::DataComponentInterpretation::component_is_scalar
                      : dealii::DataComponentInterpretation::component_is_part_of_vector);

        const std::vector<std::string> names(n_components, field.name);

        data_out.add_data_vector(dof_manager.get_dof_handler(index),
                                 solution,
                                 names,
                                 data_type);

        solution.zero_out_ghost_values();
      }

    // Build patches to linearly interpolate from higher order element degrees. Note that
    // this essentially converts the element to an equal amount of subdivisions in the
    // output. This does not make subdivisions and element degree equivalent in the
    // simulation!
    const unsigned int n_divisions = output_parameters.patch_subdivisions == 0
                                       ? degree
                                       : output_parameters.patch_subdivisions;
    data_out.build_patches(n_divisions);

    // Set some flags for data output
    dealii::DataOutBase::VtkFlags flags;
    flags.time                = sim_timer.get_time();
    flags.cycle               = sim_timer.get_increment();
    flags.print_date_and_time = true;
#ifdef PRISMS_PF_WITH_ZLIB
    flags.compression_level = dealii::DataOutBase::CompressionLevel::best_speed;
#endif
    data_out.set_flags(flags);

    // Write to file based on the user input.
    const std::string  directory = "./";
    const unsigned int increment = sim_timer.get_increment();
    const std::string &file_type = output_parameters.file_type;
    if (file_type == "vtu")
      {
        std::ostringstream increment_stream;
        increment_stream << std::setw(static_cast<int>(n_trailing_digits))
                         << std::setfill('0') << increment;
        const std::string filename =
          directory + file_prefix + "_" + increment_stream.str() + ".vtu";
        data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD);
      }
    else if (file_type == "pvtu")
      {
        data_out.write_vtu_with_pvtu_record(directory,
                                            file_prefix,
                                            increment,
                                            MPI_COMM_WORLD,
                                            n_trailing_digits);
      }
    else if (file_type == "vtk")
      {
        std::ostringstream increment_stream;
        increment_stream << std::setw(static_cast<int>(n_trailing_digits))
                         << std::setfill('0') << increment;
        const std::string filename =
          directory + file_prefix + "_" + increment_stream.str() + ".vtk";
        std::ofstream vtk_output(filename);
        data_out.write_vtk(vtk_output);
      }
    else
      {
        AssertThrow(false, UnreachableCode());
      }

    // Update the ghost values again to allow for read access
    for (unsigned int index = 0; index < field_attributes.size(); ++index)
      {
        solution_indexer.get_solution_vector(index).update_ghost_values();
      }
  }
};

PRISMS_PF_END_NAMESPACE
