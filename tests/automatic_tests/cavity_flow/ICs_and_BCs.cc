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
  [[maybe_unused]] double                   &scalar_value,
  [[maybe_unused]] double                   &vector_component_value) const
{}

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
  if (boundary_id == 3 && component == 0)
    {
      double x0 =
        0.1 * this->get_user_inputs().get_spatial_discretization().get_size()[1];

      if (point[0] <= x0)
        {
          vector_component_value = 1.0 - 0.25 *
                                           (1.0 - std::cos(M_PI * (x0 - point[0]) / x0)) *
                                           (1.0 - std::cos(M_PI * (x0 - point[0]) / x0));
        }
      else if (point[0] > x0 &&
               point[0] <
                 this->get_user_inputs().get_spatial_discretization().get_size()[1] - x0)
        {
          vector_component_value = 1.0;
        }
      else
        {
          vector_component_value =
            1.0 - 0.25 * (1.0 - std::cos(M_PI * (point[0] - (1.0 - x0)) / x0)) *
                    (1.0 - std::cos(M_PI * (point[0] - (1.0 - x0)) / x0));
        }
    }
}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE