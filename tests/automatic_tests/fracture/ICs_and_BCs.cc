// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/nonuniform_dirichlet.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/utilities.h>

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
  const double center[3] = {10.0, 12.0, 0.0};

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
      if ((std::abs(
             point[1] -
             (this->get_user_inputs().get_spatial_discretization().get_size()[1] / 2.0) +
             (0.5 * dx)) < dx) &&
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

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::set_nonuniform_dirichlet(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &boundary_id,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] number                   &scalar_value,
  [[maybe_unused]] number                   &vector_component_value) const
{
  if (index == 1)
    {
      const number time =
        this->get_user_inputs().get_temporal_discretization().get_time();

      number x = (point[0] - (vel_nom * time) - clength);
      number y =
        point[1] -
        (this->get_user_inputs().get_spatial_discretization().get_size()[1] / 2.0) +
        (dx * 0.5);
      number r      = std::sqrt((x * x) + (y * y));
      number theta  = std::atan2(y, x);
      number mu     = CIJ_base[dim][dim];
      number lambda = CIJ_base[0][0] - (2.0 * mu);
      number nu     = lambda / 2.0 / (lambda + mu);
      number kappa  = 3.0 - (4.0 * nu);

      vector_component_value = 0.5 * (KI_nom / mu) *
                               std::sqrt(0.5 * r / std::numbers::pi) *
                               (kappa - std::cos(theta));
      if (component == 0)
        {
          vector_component_value *= std::cos(0.5 * theta);
        }
      if (component == 1)
        {
          vector_component_value *= std::sin(0.5 * theta);
        }
    }
}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
