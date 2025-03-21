// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/config.h>
#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/nonuniform_dirichlet.h>

#include <cmath>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
void
customInitialCondition<dim>::set_initial_condition(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] double                   &scalar_value,
  [[maybe_unused]] double                   &vector_component_value,
  [[maybe_unused]] const userInputParameters<dim> &user_inputs) const
{
  double center[12][3] = {
    {0.1, 0.3,  0},
    {0.8, 0.7,  0},
    {0.5, 0.2,  0},
    {0.4, 0.4,  0},
    {0.3, 0.9,  0},
    {0.8, 0.1,  0},
    {0.9, 0.5,  0},
    {0.0, 0.1,  0},
    {0.1, 0.6,  0},
    {0.5, 0.6,  0},
    {1,   1,    0},
    {0.7, 0.95, 0}
  };
  double rad[12] = {12, 14, 19, 16, 11, 12, 17, 15, 20, 10, 11, 14};
  double dist    = 0.0;
  for (unsigned int i = 0; i < 12; i++)
    {
      dist = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++)
        {
          dist +=
            (point[dir] - center[i][dir] * 100) * (point[dir] - center[i][dir] * 100);
        }
      dist = std::sqrt(dist);

      scalar_value += 0.5 * (1.0 - std::tanh((dist - rad[i]) / 1.5));
    }
  scalar_value = std::min(scalar_value, 1.0);
}

template <int dim>
void
customNonuniformDirichlet<dim>::set_nonuniform_dirichlet(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &boundary_id,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] double                   &scalar_value,
  [[maybe_unused]] double                   &vector_component_value,
  [[maybe_unused]] const userInputParameters<dim> &user_inputs) const
{}

INSTANTIATE_UNI_TEMPLATE(customInitialCondition)
INSTANTIATE_UNI_TEMPLATE(customNonuniformDirichlet)

PRISMS_PF_END_NAMESPACE