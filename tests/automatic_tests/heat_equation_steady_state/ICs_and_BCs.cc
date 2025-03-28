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
  [[maybe_unused]] const unsigned int             &index,
  [[maybe_unused]] const unsigned int             &component,
  [[maybe_unused]] const dealii::Point<dim>       &point,
  [[maybe_unused]] double                         &scalar_value,
  [[maybe_unused]] double                         &vector_component_value,
  [[maybe_unused]] const userInputParameters<dim> &user_inputs) const
{
  if (index == 1)
    {
      double lx = 2.0;
      double ly = 1.0;

      scalar_value = -std::sin(M_PI * point[0] / lx) *
                     (-(M_PI / lx) * (M_PI / lx) * point[1] / ly * (1.0 - point[1] / ly) -
                      2.0 / ly / ly);
    }
}

template <int dim>
void
customNonuniformDirichlet<dim>::set_nonuniform_dirichlet(
  [[maybe_unused]] const unsigned int             &index,
  [[maybe_unused]] const unsigned int             &boundary_id,
  [[maybe_unused]] const unsigned int             &component,
  [[maybe_unused]] const dealii::Point<dim>       &point,
  [[maybe_unused]] double                         &scalar_value,
  [[maybe_unused]] double                         &vector_component_value,
  [[maybe_unused]] const userInputParameters<dim> &user_inputs) const
{}

INSTANTIATE_UNI_TEMPLATE(customInitialCondition)
INSTANTIATE_UNI_TEMPLATE(customNonuniformDirichlet)

PRISMS_PF_END_NAMESPACE