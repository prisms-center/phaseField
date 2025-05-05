// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/nonuniform_dirichlet.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <cmath>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
void
customInitialCondition<dim>::set_initial_condition(
  [[maybe_unused]] const unsigned int             &index,
  [[maybe_unused]] const unsigned int             &component,
  [[maybe_unused]] const dealii::Point<dim>       &point,
  [[maybe_unused]] double                         &scalar_value,
  [[maybe_unused]] double                         &vector_component_value,
  [[maybe_unused]] const userInputParameters<dim> &user_inputs) const
{
  const double center[3] = {10.0, 12.0, 0.0};
  const double dx        = user_inputs.spatial_discretization.size[0] /
                    ((double) user_inputs.spatial_discretization.subdivisions[0] /
                     std::pow(2.0, user_inputs.spatial_discretization.global_refinement));
  const double clength =
    user_inputs.user_constants.get_model_constant_double("cracklength");

  double dist = 0.0;
  for (unsigned int dir = 0; dir < dim; dir++)
    {
      dist += (point[dir] - center[dir]) * (point[dir] - center[dir]);
    }
  dist = std::sqrt(dist);

  scalar_value           = 0.0;
  vector_component_value = 0.0;

  if (index == 0)
    {
      if ((std::abs(point[1] - user_inputs.spatial_discretization.size[1] / 2.0 +
                    0.5 * dx) < dx) &&
          (point[0] < clength))
        {
          scalar_value = 1.0;
        }
    }
  if (index == 3 || index == 4)
    {
      scalar_value = 1.0;
    }
}

template <unsigned int dim, typename number>
void
customNonuniformDirichlet<dim, number>::set_nonuniform_dirichlet(
  [[maybe_unused]] const unsigned int             &index,
  [[maybe_unused]] const unsigned int             &boundary_id,
  [[maybe_unused]] const unsigned int             &component,
  [[maybe_unused]] const dealii::Point<dim>       &point,
  [[maybe_unused]] number                         &scalar_value,
  [[maybe_unused]] number                         &vector_component_value,
  [[maybe_unused]] const userInputParameters<dim> &user_inputs) const
{}

template class customInitialCondition<1>;
template class customInitialCondition<2>;
template class customInitialCondition<3>;

template class customNonuniformDirichlet<1, double>;
template class customNonuniformDirichlet<2, double>;
template class customNonuniformDirichlet<3, double>;
template class customNonuniformDirichlet<1, float>;
template class customNonuniformDirichlet<2, float>;
template class customNonuniformDirichlet<3, float>;

PRISMS_PF_END_NAMESPACE
