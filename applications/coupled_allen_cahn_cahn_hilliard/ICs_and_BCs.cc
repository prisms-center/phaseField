// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/nonuniform_dirichlet.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include <algorithm>
#include <cmath>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::set_initial_condition(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] number                   &scalar_value,
  [[maybe_unused]] number                   &vector_component_value) const
{
  double dist;

  if (index == 0)
    {
      scalar_value = matrix_concentration;
    }

  dist = 0.0;
  for (unsigned int dir = 0; dir < dim; dir++)
    {
      dist += (point[dir] - center1[dir]) * (point[dir] - center1[dir]);
    }
  dist = std::sqrt(dist);

  // Initial condition for the concentration field
  if (index == 0)
    {
      scalar_value += 0.5 * (0.125) * (1.0 - std::tanh((dist - radius1) / (1.0)));
    }
  else
    {
      scalar_value += 0.5 * (1.0 - std::tanh((dist - radius1) / (1.0)));
    }

  dist = 0.0;
  for (unsigned int dir = 0; dir < dim; dir++)
    {
      dist += (point[dir] - center2[dir]) * (point[dir] - center2[dir]);
    }
  dist = std::sqrt(dist);

  // Initial condition for the concentration field
  if (index == 0)
    {
      scalar_value += 0.5 * (0.125) * (1.0 - std::tanh((dist - radius2) / (1.0)));
    }
  else
    {
      scalar_value += 0.5 * (1.0 - std::tanh((dist - radius2) / (1.0)));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::set_nonuniform_dirichlet(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &boundary_id,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] number                   &scalar_value,
  [[maybe_unused]] number                   &vector_component_value) const
{}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
