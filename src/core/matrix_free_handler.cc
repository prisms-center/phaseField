// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/quadrature.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/timer.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <memory>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, typename number>
MatrixFreeHandler<dim, number>::MatrixFreeHandler()
  : matrix_free_object(
      std::make_shared<
        dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>>())
{
  additional_data.tasks_parallel_scheme =
    dealii::MatrixFree<dim,
                       number>::AdditionalData::TasksParallelScheme::partition_partition;

  // TODO (landinjm): This should be done according to the residual flags to prevent
  // excess data being evaluated for the shape functions. Note that this applies to all
  // PDEs, so it might now be too wasteful.
  additional_data.mapping_update_flags =
    (dealii::update_values | dealii::update_gradients | dealii::update_JxW_values |
     dealii::update_quadrature_points);
}

template <unsigned int dim, typename number>
void
MatrixFreeHandler<dim, number>::reinit(
  const dealii::Mapping<dim>              &mapping,
  const dealii::DoFHandler<dim>           &dof_handler,
  const dealii::AffineConstraints<number> &constraint,
  const dealii::Quadrature<1>             &quad)
{
  matrix_free_object->reinit(mapping, dof_handler, constraint, quad, additional_data);
  is_initialized = true;
}

template <unsigned int dim, typename number>
void
MatrixFreeHandler<dim, number>::reinit(
  const dealii::Mapping<dim>                                   &mapping,
  const std::vector<const dealii::DoFHandler<dim> *>           &dof_handler,
  const std::vector<const dealii::AffineConstraints<number> *> &constraint,
  const dealii::Quadrature<1>                                  &quad)
{
  matrix_free_object->reinit(mapping, dof_handler, constraint, quad, additional_data);
  is_initialized = true;
}

template <unsigned int dim, typename number>
void
MatrixFreeHandler<dim, number>::reinit(
  const dealii::Mapping<dim>                                   &mapping,
  const std::vector<const dealii::DoFHandler<dim> *>           &dof_handler,
  const std::vector<const dealii::AffineConstraints<number> *> &constraint,
  const std::vector<dealii::Quadrature<1>>                     &quad)
{
  matrix_free_object->reinit(mapping, dof_handler, constraint, quad, additional_data);
  is_initialized = true;
}

template <unsigned int dim, typename number>
std::shared_ptr<dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>>
MatrixFreeHandler<dim, number>::get_matrix_free() const
{
  Assert(is_initialized, dealii::ExcNotInitialized());
  return matrix_free_object;
}

template <unsigned int dim, typename number>
MatrixFreeContainer<dim, number>::MatrixFreeContainer(MGInfo<dim> &mg_info)
  : matrix_free()
  , multigrid_matrix_free(0, 0)
{
  // If we have multigrid, we need to set the min and max levels, and reinitialize the
  // multigrid matrix-free handler.
  if (mg_info.has_multigrid())
    {
      min_level = mg_info.get_mg_min_level();
      max_level = mg_info.get_mg_max_level();
      multigrid_matrix_free.resize(min_level, max_level);
    }
}

template <unsigned int dim, typename number>
template <unsigned int degree, unsigned int quad_dim>
void
MatrixFreeContainer<dim, number>::reinit(
  const dealii::Mapping<dim>                   &mapping,
  const DofHandler<dim>                        &dof_container,
  const ConstraintHandler<dim, degree, number> &constraint_container,
  const dealii::Quadrature<quad_dim>           &quad)
{
  Timer::start_section("reinitialize matrix-free objects");

  ConditionalOStreams::pout_base() << "initializing matrix-free object...\n"
                                   << std::flush;
  // Reinit the matrix-free object
  matrix_free.reinit(mapping,
                     dof_container.get_dof_handlers(),
                     constraint_container.get_constraints(),
                     quad);

  // Reinit the multigrid matrix-free objects if we have multigrid
  if (multigrid_matrix_free.n_levels() > 1)
    {
      for (unsigned int level = min_level; level <= max_level; ++level)
        {
          ConditionalOStreams::pout_base()
            << "initializing multgrid matrix-free object at level " << level << "...\n"
            << std::flush;
          multigrid_matrix_free[level].reinit(mapping,
                                              dof_container.get_mg_dof_handlers(level),
                                              constraint_container.get_mg_constraints(
                                                level),
                                              quad);
        }
    }
  Timer::end_section("reinitialize matrix-free objects");
}

template <unsigned int dim, typename number>
std::shared_ptr<dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>>
MatrixFreeContainer<dim, number>::get_matrix_free() const
{
  return matrix_free.get_matrix_free();
}

template <unsigned int dim, typename number>
std::shared_ptr<dealii::MatrixFree<dim, float, dealii::VectorizedArray<float>>>
MatrixFreeContainer<dim, number>::get_mg_matrix_free(unsigned int level) const
{
  Assert(multigrid_matrix_free.n_levels() > level,
         dealii::ExcMessage("The multigrid matrix-free object does not contain level = " +
                            std::to_string(level)));
  return multigrid_matrix_free[level].get_matrix_free();
}

#include "core/matrix_free_handler.inst"

PRISMS_PF_END_NAMESPACE
