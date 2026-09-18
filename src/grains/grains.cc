// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/grains/grains.h>

PRISMS_PF_BEGIN_NAMESPACE

// ============================================================================
// Methods for SimplifiedGrainRepresentation
// ============================================================================

template <unsigned int dim>
SimplifiedGrainRepresentation<dim>::SimplifiedGrainRepresentation(
  const GrainSet<dim> &grain_set)
  : grain_id(grain_set.get_grain_index())
  , order_parameter_id(grain_set.get_order_parameter_index())
  , old_order_parameter_id(order_parameter_id)
  , distance_to_neighbor_sharing_op(0.0)
{
  // Calculate the centroid assuming that the elements are rectangular and with
  // no weighting based on the actual value of the field
  std::vector<std::vector<dealii::Point<dim>>> vertex_list = grain_set.get_vertex_list();

  double                 grain_volume = 0.0;
  dealii::Tensor<1, dim> centroid;

  for (auto &vertex : vertex_list)
    {
      double             cell_volume = 1.0;
      dealii::Point<dim> cell_center;

      unsigned int opposite_corner_index = 0;
      if (dim == 2)
        {
          opposite_corner_index = 3;
        }
      else
        {
          opposite_corner_index = 7;
        }

      for (unsigned int dimension = 0; dimension < dim; dimension++)
        {
          cell_volume *=
            (vertex[opposite_corner_index][dimension] - vertex[0][dimension]);
          cell_center(dimension) =
            (vertex[opposite_corner_index][dimension] + vertex[0][dimension]) / 2.0;
        }

      for (unsigned int dimension = 0; dimension < dim; dimension++)
        {
          centroid[dimension] += cell_volume * cell_center(dimension);
        }

      grain_volume += cell_volume;
    }

  centroid /= grain_volume;

  for (unsigned int dimension = 0; dimension < dim; dimension++)
    {
      center(dimension) = centroid[dimension];
    }

  // Calculate the radius as the largest distance from the centroid to one of
  // the vertices
  radius = 0.0;
  for (auto &vertex : vertex_list)
    {
      for (unsigned int vertex_index = 0;
           vertex_index < dealii::Utilities::fixed_power<dim>(2.0);
           vertex_index++)
        {
          if (vertex[vertex_index].distance(center) > radius)
            {
              radius = vertex[vertex_index].distance(center);
            }
        }
    }
}

template <unsigned int dim>
dealii::Point<dim>
SimplifiedGrainRepresentation<dim>::get_center() const
{
  return center;
}

template <unsigned int dim>
double
SimplifiedGrainRepresentation<dim>::get_radius() const
{
  return radius;
}

template <unsigned int dim>
unsigned int
SimplifiedGrainRepresentation<dim>::get_grain_id() const
{
  return grain_id;
}

template <unsigned int dim>
void
SimplifiedGrainRepresentation<dim>::set_grain_id(unsigned int _grain_id)
{
  grain_id = _grain_id;
}

template <unsigned int dim>
unsigned int
SimplifiedGrainRepresentation<dim>::get_order_parameter_id() const
{
  return order_parameter_id;
}

template <unsigned int dim>
void
SimplifiedGrainRepresentation<dim>::set_order_parameter_id(
  unsigned int _order_parameter_id)
{
  order_parameter_id = _order_parameter_id;
}

template <unsigned int dim>
unsigned int
SimplifiedGrainRepresentation<dim>::get_old_order_parameter_id() const
{
  return old_order_parameter_id;
}

template <unsigned int dim>
void
SimplifiedGrainRepresentation<dim>::set_distance_to_neighbor(double dist)
{
  distance_to_neighbor_sharing_op = dist;
}

template <unsigned int dim>
double
SimplifiedGrainRepresentation<dim>::get_distance_to_neighbor() const
{
  return distance_to_neighbor_sharing_op;
}

#include "grains/grains.inst"

PRISMS_PF_END_NAMESPACE
