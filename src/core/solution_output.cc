#include <prismspf/core/solution_output.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

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
  // TODO (landinjm): Make this a user input parameter so they can select between
  // compression levels. Best speed can be the default
  flags.compression_level = dealii::DataOutBase::CompressionLevel::best_speed;
#endif
  data_out.set_flags(flags);

  // Write to file based on the user input.
  // TODO (landinjm): actually write stuff according to user input.
  data_out.write_vtu_with_pvtu_record("./",
                                      name,
                                      user_inputs.temporal_discretization.increment,
                                      MPI_COMM_WORLD,
                                      n_trailing_digits);
}

template <int dim, typename number>
solutionOutput<dim, number>::solutionOutput(
  const std::map<unsigned int, VectorType *>         &solution_set,
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
  for (const auto &[index, variable] : *user_inputs.var_attributes)
    {
      auto *solution = solution_set.at(index);
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

  // Write to file based on the user input.
  // TODO (landinjm): actually write stuff according to user input.
  data_out.write_vtu_with_pvtu_record("./",
                                      name,
                                      user_inputs.temporal_discretization.increment,
                                      MPI_COMM_WORLD,
                                      n_trailing_digits);
}

template class solutionOutput<1, float>;
template class solutionOutput<1, double>;
template class solutionOutput<2, float>;
template class solutionOutput<2, double>;
template class solutionOutput<3, float>;
template class solutionOutput<3, double>;

PRISMS_PF_END_NAMESPACE