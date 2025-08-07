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
  if constexpr (dim == 2)
    {
      const double x = point[0];
      const double y = point[1];

      if (index == 0 || index == 1)
        {
          if (component == 0)
            {
              vector_component_value = std::sin(x) * std::cos(y);
            }
          if (component == 1)
            {
              vector_component_value = -std::cos(x) * std::sin(y);
            }
        }
    }
  else if constexpr (dim == 3)
    {
      const double x = point[0];
      const double y = point[1];
      const double z = point[2];

      if (index == 0 || index == 1)
        {
          if (component == 0)
            {
              vector_component_value = std::sin(x) * std::cos(y) * std::cos(z);
            }
          if (component == 1)
            {
              vector_component_value = -std::cos(x) * std::sin(y) * std::cos(z);
            }
          if (component == 2)
            {
              vector_component_value = 0;
            }
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
