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
customPDE<dim, degree, number>::set_initial_condition(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] double                   &scalar_value,
  [[maybe_unused]] double                   &vector_component_value) const
{
  const double center[3] = {10.0, 12.0, 0.0};
  const double dx =
    this->get_user_inputs().spatial_discretization.size[0] /
    double(this->get_user_inputs().spatial_discretization.subdivisions[0]) /
    std::pow(2.0, this->get_user_inputs().spatial_discretization.global_refinement);
  const double clength =
    this->get_user_inputs().user_constants.get_model_constant_double("cracklength");

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
      if ((std::abs(point[1] -
                    (this->get_user_inputs().spatial_discretization.size[1] / 2.0) +
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
customPDE<dim, degree, number>::set_nonuniform_dirichlet(
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
        this->get_user_inputs().temporal_discretization.get_current_time();

      const number dx =
        this->get_user_inputs().spatial_discretization.size[0] /
        number(this->get_user_inputs().spatial_discretization.subdivisions[0]) /
        std::pow(2.0, this->get_user_inputs().spatial_discretization.global_refinement);

      const number clength =
        this->get_user_inputs().user_constants.get_model_constant_double("cracklength");
      const dealii::Tensor<2, voigt_tensor_size<dim>, number> CIJ_base =
        this->get_user_inputs().user_constants.get_model_constant_elasticity_tensor(
          "CIJ_base");
      const number KI_nom =
        this->get_user_inputs().user_constants.get_model_constant_double("KI_nom");
      const number vel_nom =
        this->get_user_inputs().user_constants.get_model_constant_double("vel_nom");

      number x = (point[0] - (vel_nom * time) - clength);
      number y = point[1] -
                 (this->get_user_inputs().spatial_discretization.size[1] / 2.0) +
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

INSTANTIATE_TRI_TEMPLATE(customPDE)

PRISMS_PF_END_NAMESPACE
