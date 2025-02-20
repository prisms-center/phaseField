// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <cmath>
#include <config.h>
#include <core/boundary_conditions/nonuniform_dirichlet.h>

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
      if (boundary_id == 0 || boundary_id == 1)
        {
          scalar_value = std::sin(point[1] * M_PI);
        }
    }
}

INSTANTIATE_UNI_TEMPLATE(customNonuniformDirichlet)