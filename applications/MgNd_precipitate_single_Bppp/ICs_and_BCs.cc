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
  double center_x = 0.5 * get_user_inputs().get_spatial_discretization().get_size()[0];
  double center_y = 0.5 * get_user_inputs().get_spatial_discretization().get_size()[1];
  double center_z = 0.5 * get_user_inputs().get_spatial_discretization().get_size()[2];
  dealii::Point<3> center(center_x, center_y, center_z);

  // Initial condition parameters
  double x_denom = (1.0) * (1.0);
  double y_denom = (8.0) * (8.0);
  double z_denom = (8.0) * (8.0);

  double initial_interface_coeff = 0.08;
  double initial_radius          = 1.0;
  double c_matrix                = 1.0e-6;
  double c_precip                = 0.14;

  // set result equal to the structural order parameter initial condition
  number              r = 0.0;
  std::vector<double> ellipsoid_denoms;
  ellipsoid_denoms.push_back(x_denom);
  ellipsoid_denoms.push_back(y_denom);
  ellipsoid_denoms.push_back(z_denom);

  for (unsigned int i = 0; i < dim; i++)
    {
      r += (point(i) - center[i]) * (point(i) - center[i]) / ellipsoid_denoms[i];
    }
  r = std::sqrt(r);

  if (index == 0)
    {
      scalar_value =
        0.5 * (c_precip - c_matrix) *
          (1.0 - std::tanh((r - initial_radius) / (initial_interface_coeff))) +
        c_matrix;
    }
  else if (index == 1)
    {
      scalar_value = 0.0;
    }
  else if (index == 2)
    {
      scalar_value =
        0.5 * (1.0 - std::tanh((r - initial_radius) / (initial_interface_coeff)));
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
