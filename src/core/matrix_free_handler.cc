// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/matrix_free_handler.h>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim, typename number>
matrixfreeHandler<dim, number>::matrixfreeHandler(
  const userInputParameters<dim> &_user_inputs)
  : user_inputs(_user_inputs)
  , matrix_free_object(std::make_shared<dealii::MatrixFree<dim, number>>())
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

template <int dim, typename number>
void
matrixfreeHandler<dim, number>::reinit(
  const dealii::Mapping<dim>              &mapping,
  const dealii::DoFHandler<dim>           &dof_handler,
  const dealii::AffineConstraints<number> &constraint,
  const dealii::Quadrature<1>             &quad)
{
  matrix_free_object->reinit(mapping, dof_handler, constraint, quad, additional_data);
}

template <int dim, typename number>
void
matrixfreeHandler<dim, number>::reinit(
  const dealii::Mapping<dim>                                   &mapping,
  const std::vector<const dealii::DoFHandler<dim> *>           &dof_handler,
  const std::vector<const dealii::AffineConstraints<number> *> &constraint,
  const dealii::Quadrature<1>                                  &quad)
{
  matrix_free_object->reinit(mapping, dof_handler, constraint, quad, additional_data);
}

template <int dim, typename number>
void
matrixfreeHandler<dim, number>::reinit(
  const dealii::Mapping<dim>                                   &mapping,
  const std::vector<const dealii::DoFHandler<dim> *>           &dof_handler,
  const std::vector<const dealii::AffineConstraints<number> *> &constraint,
  const std::vector<dealii::Quadrature<1>>                     &quad)
{
  matrix_free_object->reinit(mapping, dof_handler, constraint, quad, additional_data);
}

template <int dim, typename number>
std::shared_ptr<dealii::MatrixFree<dim, number>>
matrixfreeHandler<dim, number>::get_matrix_free() const
{
  return matrix_free_object;
}

template class matrixfreeHandler<1, double>;
template class matrixfreeHandler<2, double>;
template class matrixfreeHandler<3, double>;
template class matrixfreeHandler<1, float>;
template class matrixfreeHandler<2, float>;
template class matrixfreeHandler<3, float>;

PRISMS_PF_END_NAMESPACE