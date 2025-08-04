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
  double center[1][3] = {
    {0.0, 0.0, 0.0}
  };
  double rad[1] = {5};
  double dist;

  // Initial condition for the concentration field
  if (index == 0)
    {
      scalar_value = U0;
    }
  // Initial condition for the order parameter field
  else if (index == 1)
    {
      // Initial condition for the order parameter field
      for (unsigned int i = 0; i < 1; i++)
        {
          dist = 0.0;
          for (unsigned int dir = 0; dir < dim; dir++)
            {
              dist +=
                (point[dir] -
                 center[i][dir] *
                   this->get_user_inputs().get_spatial_discretization().get_size()[dir]) *
                (point[dir] -
                 center[i][dir] *
                   this->get_user_inputs().get_spatial_discretization().get_size()[dir]);
            }
          dist = std::sqrt(dist);

          scalar_value += (-std::tanh((dist - rad[i]) / (sqrt(2))));
        }
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