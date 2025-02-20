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
  [[maybe_unused]] const uint               &index,
  [[maybe_unused]] const uint               &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] double                   &scalar_value,
  [[maybe_unused]] double                   &vector_component_value) const
{}

template <int dim>
void
customNonuniformDirichlet<dim>::set_nonuniform_dirichlet(
  [[maybe_unused]] const uint               &index,
  [[maybe_unused]] const uint               &boundary_id,
  [[maybe_unused]] const uint               &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] double                   &scalar_value,
  [[maybe_unused]] double                   &vector_component_value) const
{
  if (index == 0)
    {
      scalar_value = std::sin(2.0 * M_PI * (point[0] + point[1]));
    }
}

INSTANTIATE_UNI_TEMPLATE(customInitialCondition)
INSTANTIATE_UNI_TEMPLATE(customNonuniformDirichlet)

PRISMS_PF_END_NAMESPACE