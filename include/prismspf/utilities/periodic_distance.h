// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/point.h>

#include <prismspf/user_inputs/spatial_discretization.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/utilities.h>

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
         const RectangularMesh<dim>     &rectangular_mesh)
{
  if (rectangular_mesh.periodic_directions.empty())
    {
      return point1.distance(point2);
    }
  // else
  real dist = real(0.0);
  for (unsigned int d = 0; d < dim; ++d)
    {
      using std::min;
      real delta = point2[d] - point1[d];
      if (rectangular_mesh.periodic_directions.contains(d))
        {
          const real length      = rectangular_mesh.size[d];
          const real half_length = length / 2.0;
          delta                  = pmod(delta - half_length, length) - half_length;
        }
      dist += delta * delta;
    }
  using std::sqrt;
  return sqrt(dist);
}

PRISMS_PF_END_NAMESPACE
