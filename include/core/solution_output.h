#ifndef solution_output_h
#define solution_output_h

#include <deal.II/base/mpi.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/numerics/data_out.h>

#include <string>

/**
 * \brief Class that outputs a passed solution to vtu, vtk, or pvtu
 */
template <int dim>
class solutionOutput
{
public:
  /**
   * \brief Constructor for a single field that must be output.
   */
  solutionOutput(const dealii::LinearAlgebra::distributed::Vector<double> &solution,
                 const dealii::DoFHandler<dim>                            &dof_handler,
                 const uint                                               &degree,
                 const std::string                                        &name,
                 const uint                                               &increment);

private:
};

template <int dim>
solutionOutput<dim>::solutionOutput(
  const dealii::LinearAlgebra::distributed::Vector<double> &solution,
  const dealii::DoFHandler<dim>                            &dof_handler,
  const uint                                               &degree,
  const std::string                                        &name,
  const uint                                               &increment)
{
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
  flags.compression_level = dealii::DataOutBase::CompressionLevel::best_speed;
  data_out.set_flags(flags);
  data_out.write_vtu_with_pvtu_record("./", name, increment, MPI_COMM_WORLD, 3);
}

#endif