// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/point.h>

#include <prismspf/user_inputs/boundary_parameters.h>
#include <prismspf/user_inputs/spatial_discretization.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Calculate the distance between two points considering periodic boundaries.
 * Note: templates must be provided explicitly because dealii::Point uses a signed int for
 * dimension, disallowing automatic deduction.
 */
template <unsigned int dim, typename real>
real
distance(const dealii::Point<dim, real> &point1,
         const dealii::Point<dim, real> &point2,
         const UserInputParameters<dim> &user_inputs)
{
  const BoundaryParameters<dim> &boundary_params = user_inputs.get_boundary_parameters();
  const SpatialDiscretization<dim> &spatial_discretization =
    user_inputs.get_spatial_discretization();
  if (!boundary_params.has_periodic_boundaries())
    {
      return point1.distance(point2);
    }
  // else
  real dist = real(0.0);
  for (unsigned int d = 0; d < dim; ++d)
    {
      using std::min;
      real delta = point2[d] - point1[d];
      if (boundary_params.get_periodicity()[d])
        {
          using std::fmod;
          const real half_box_length = spatial_discretization.get_size()[d] / 2.0;
          delta = fmod(delta - half_box_length, spatial_discretization.get_size()[d]) -
                  half_box_length;
        }
      dist += delta * delta;
    }
  using std::sqrt;
  return sqrt(dist);
}

PRISMS_PF_END_NAMESPACE
