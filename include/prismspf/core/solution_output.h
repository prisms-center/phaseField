// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef solution_output_h
#define solution_output_h

#include <deal.II/base/mpi.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/numerics/data_out.h>

#include <prismspf/config.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Class that outputs a passed solution to vtu, vtk, or pvtu
 */
template <int dim, typename number = double>
class solutionOutput
{
public:
  using VectorType = dealii::LinearAlgebra::distributed::Vector<number>;

  /**
   * \brief Constructor for a single field that must be output.
   */
  solutionOutput(const VectorType               &solution,
                 const dealii::DoFHandler<dim>  &dof_handler,
                 const unsigned int             &degree,
                 const std::string              &name,
                 const userInputParameters<dim> &user_inputs);

  /**
   * \brief Constructor for a multiple fields that must be output.
   */
  solutionOutput(const std::unordered_map<std::pair<unsigned int, dependencyType>,
                                          VectorType *,
                                          pairHash>                 &solution_set,
                 const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers,
                 const unsigned int                                 &degree,
                 const std::string                                  &name,
                 const userInputParameters<dim>                     &user_inputs);

private:
};

template <int dim, typename number>
solutionOutput<dim, number>::solutionOutput(const VectorType               &solution,
                                            const dealii::DoFHandler<dim>  &dof_handler,
                                            const unsigned int             &degree,
                                            const std::string              &name,
                                            const userInputParameters<dim> &user_inputs)
{
  // Some stuff to determine the actual name of the output file.
  const auto n_trailing_digits = static_cast<unsigned int>(
    std::floor(std::log10(user_inputs.temporal_discretization.total_increments)) + 1);

  // Init data out
  dealii::DataOut<dim> data_out;

  // Add data vector
  data_out.add_data_vector(dof_handler, solution, name);

  // Build patches to linearly interpolate from higher order element degrees. Note that
  // this essentially converts the element to an equal amount of subdivisions in the
  // output. This does not make subdivisions and element degree equivalent in the
  // simulation!
  data_out.build_patches(degree);

  // Set some flags for data output
  dealii::DataOutBase::VtkFlags flags;
  flags.time                = user_inputs.temporal_discretization.time;
  flags.cycle               = user_inputs.temporal_discretization.increment;
  flags.print_date_and_time = true;
#ifdef PRISMS_PF_WITH_ZLIB
  flags.compression_level = dealii::DataOutBase::CompressionLevel::best_speed;
#endif
  data_out.set_flags(flags);

  // Write to file based on the user input. TODO: actually write stuff according to user
  // input.
  data_out.write_vtu_with_pvtu_record("./",
                                      name,
                                      user_inputs.temporal_discretization.increment,
                                      MPI_COMM_WORLD,
                                      n_trailing_digits);
}

template <int dim, typename number>
solutionOutput<dim, number>::solutionOutput(
  const std::unordered_map<std::pair<unsigned int, dependencyType>,
                           VectorType *,
                           pairHash>                 &solution_set,
  const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers,
  const unsigned int                                 &degree,
  const std::string                                  &name,
  const userInputParameters<dim>                     &user_inputs)
{
  // Some stuff to determine the actual name of the output file.
  const auto n_trailing_digits = static_cast<unsigned int>(
    std::floor(std::log10(user_inputs.temporal_discretization.total_increments)) + 1);

  // Init data out
  dealii::DataOut<dim> data_out;

  // Add data vectors
  for (const auto &[index, variable] : user_inputs.var_attributes)
    {
      auto *solution = solution_set.at(std::make_pair(index, dependencyType::NORMAL));
      solution->update_ghost_values();

      // Mark field as SCALAR/VECTOR
      const bool         is_scalar    = variable.field_type == fieldType::SCALAR;
      const unsigned int n_components = is_scalar ? 1 : dim;

      const std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
        dataType(n_components,
                 is_scalar
                   ? dealii::DataComponentInterpretation::component_is_scalar
                   : dealii::DataComponentInterpretation::component_is_part_of_vector);

      const std::vector<std::string> names(n_components, variable.name.c_str());

      data_out.add_data_vector(*(dof_handlers.at(index)), *solution, names, dataType);

      solution->zero_out_ghost_values();
    }

  // Build patches to linearly interpolate from higher order element degrees. Note that
  // this essentially converts the element to an equal amount of subdivisions in the
  // output. This does not make subdivisions and element degree equivalent in the
  // simulation!
  data_out.build_patches(degree);

  // Set some flags for data output
  dealii::DataOutBase::VtkFlags flags;
  flags.time                = user_inputs.temporal_discretization.time;
  flags.cycle               = user_inputs.temporal_discretization.increment;
  flags.print_date_and_time = true;
#ifdef PRISMS_PF_WITH_ZLIB
  flags.compression_level = dealii::DataOutBase::CompressionLevel::best_speed;
#endif
  data_out.set_flags(flags);

  // Write to file based on the user input. TODO: actually write stuff according to user
  // input.
  data_out.write_vtu_with_pvtu_record("./",
                                      name,
                                      user_inputs.temporal_discretization.increment,
                                      MPI_COMM_WORLD,
                                      n_trailing_digits);
}

PRISMS_PF_END_NAMESPACE

#endif