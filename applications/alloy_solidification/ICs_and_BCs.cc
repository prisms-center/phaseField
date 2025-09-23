// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/nonuniform_dirichlet.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

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
  const std::array<double, 3> center = {
    {0.0, 0.0, 0.0}
  };
  const double rad  = 5.0;
  double       dist = 0.0;

  if (index == 0)
    {
      // For the concentration field, we just set the initial condition
      // to have uniform undercooling everywhere.
      scalar_value = U0;
    }
  else if (index == 1)
    {
      // For the order parameter, we just place a small seed. Note that
      // the order parameter ranges from -1 to 1 in thi model.
      for (unsigned int dir = 0; dir < dim; dir++)
        {
          dist += (point[dir] -
                   center[dir] *
                     get_user_inputs().get_spatial_discretization().get_size()[dir]) *
                  (point[dir] -
                   center[dir] *
                     get_user_inputs().get_spatial_discretization().get_size()[dir]);
        }
      dist = std::sqrt(dist);

      scalar_value = -std::tanh((dist - rad) / std::sqrt(2));
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
