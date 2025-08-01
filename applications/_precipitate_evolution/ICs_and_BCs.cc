// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/nonuniform_dirichlet.h>

#include <prismspf/config.h>

#include <cmath>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
void
customInitialCondition<dim>::set_initial_condition(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] number                   &scalar_value,
  [[maybe_unused]] number                   &vector_component_value) const
{
  double center[4][3] = {
    {1.0 / 3.0, 1.0 / 3.0, 0.5},
    {2.0 / 3.0, 2.0 / 3.0, 0.5},
    {3.0 / 4.0, 1.0 / 4.0, 0.5},
    {1.0 / 4.0, 3.0 / 4,   0.5}
  };
  double rad[4]         = {40.0 / 16.0, 40.0 / 16.0, 40.0 / 16.0, 40.0 / 16.0};
  double orientation[4] = {1, 1, 2, 3};
  double dx             = 40.0 / (3.0 / std::pow(2.0, 5));
  double dist           = 0.0;

  if (index == 0)
    {
      scalar_value = 0.04;
    }

  for (unsigned int i = 0; i < 4; i++)
    {
      dist = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++)
        {
          dist +=
            (point[dir] - center[i][dir] * 40.0) * (point[dir] - center[i][dir] * 40.0);
        }
      dist = std::sqrt(dist);

      if (index == orientation[i])
        {
          scalar_value += 0.5 * (1.0 - std::tanh((dist - rad[i]) / (dx)));
        }
    }

  if (index == 4)
    {
      vector_component_value = 0.0;
    }
}

template <int dim>
void
customNonuniformDirichlet<dim>::set_nonuniform_dirichlet(
  [[maybe_unused]] const unsigned int             &index,
  [[maybe_unused]] const unsigned int             &boundary_id,
  [[maybe_unused]] const unsigned int             &component,
  [[maybe_unused]] const dealii::Point<dim>       &point,
  [[maybe_unused]] number                         &scalar_value,
  [[maybe_unused]] number                         &vector_component_value,
  [[maybe_unused]] const UserInputParameters<dim> &user_inputs) const
{}

INSTANTIATE_UNI_TEMPLATE(customInitialCondition)
INSTANTIATE_UNI_TEMPLATE(customNonuniformDirichlet)

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
