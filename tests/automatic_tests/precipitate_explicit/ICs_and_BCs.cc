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
{
  const double x_length =
    this->get_user_inputs().get_spatial_discretization().get_size()[0];
  const double refinement =
    this->get_user_inputs().get_spatial_discretization().get_global_refinement();
  const double radius_size_factor = 16.0;
  const double concentration      = 0.04;

  // Oreintation of the four precipitates
  std::vector<unsigned int> orientation = {1, 1, 2, 3};

  // Centers of the precipitates in fractional coordinates
  std::vector<std::vector<double>> center = {
    {1.0 / 3.0, 1.0 / 3.0, 0.5},
    {2.0 / 3.0, 2.0 / 3.0, 0.5},
    {3.0 / 4.0, 1.0 / 4.0, 0.5},
    {1.0 / 4.0, 3.0 / 4.0, 0.5}
  };
  std::vector<double> rad = {x_length / radius_size_factor,
                             x_length / radius_size_factor,
                             x_length / radius_size_factor,
                             x_length / radius_size_factor};

  double dx = x_length / (3.0 * std::pow(2.0, refinement));

  if (index == 0)
    {
      scalar_value = concentration;
    }
  if (index == 1 || index == 2 || index == 3)
    {
      scalar_value = 0.0;
      for (unsigned int i = 0; i < orientation.size(); i++)
        {
          double dist = 0.0;
          for (unsigned int dir = 0; dir < dim; dir++)
            {
              dist += (point[dir] - center[i].at(dir) * x_length) *
                      (point[dir] - center[i].at(dir) * x_length);
            }
          dist = std::sqrt(dist);

          if (index == orientation[i])
            {
              scalar_value += 0.5 * (1.0 - std::tanh((dist - rad[i]) / dx));
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

INSTANTIATE_TRI_TEMPLATE(CustomPDE)

PRISMS_PF_END_NAMESPACE
