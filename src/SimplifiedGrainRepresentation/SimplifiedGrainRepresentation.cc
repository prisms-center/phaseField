#include "../../include/SimplifiedGrainRepresentation.h"

// ============================================================================
// Methods for SimplifiedGrainRepresentation
// ============================================================================

template <int dim>
SimplifiedGrainRepresentation<dim>::SimplifiedGrainRepresentation(
  const GrainSet<dim> &grain_set)
{
  grain_id                        = grain_set.getGrainIndex();
  order_parameter_id              = grain_set.getOrderParameterIndex();
  old_order_parameter_id          = order_parameter_id;
  distance_to_neighbor_sharing_op = 0.0;

  // Calculate the centroid assuming that the elements are rectangular and with
  // no weighting based on the actual value of the field
  std::vector<std::vector<dealii::Point<dim>>> vertex_list = grain_set.getVertexList();

  double                 grain_volume = 0.0;
  dealii::Tensor<1, dim> centroid;

  for (unsigned int c = 0; c < vertex_list.size(); c++)
    {
      double             cell_volume = 1.0;
      dealii::Point<dim> cell_center;

      unsigned int opposite_corner_index;
      if (dim == 2)
        {
          opposite_corner_index = 3;
        }
      else
        {
          opposite_corner_index = 7;
        }

      for (unsigned int d = 0; d < dim; d++)
        {
          cell_volume *=
            (vertex_list.at(c)[opposite_corner_index][d] - vertex_list.at(c)[0][d]);
          cell_center(d) =
            (vertex_list.at(c)[opposite_corner_index][d] + vertex_list.at(c)[0][d]) / 2.0;
        }

      for (unsigned int d = 0; d < dim; d++)
        {
          centroid[d] += cell_volume * cell_center(d);
        }

      grain_volume += cell_volume;
    }

  centroid /= grain_volume;

  for (unsigned int d = 0; d < dim; d++)
    {
      center(d) = centroid[d];
    }

  // Calculate the radius as the largest distance from the centroid to one of
  // the vertices
  radius = 0.0;
  for (unsigned int c = 0; c < vertex_list.size(); c++)
    {
      for (unsigned int v = 0; v < dealii::Utilities::fixed_power<dim>(2.0); v++)
        {
          if (vertex_list[c][v].distance(center) > radius)
            {
              radius = vertex_list[c][v].distance(center);
            }
        }
    }
}

template <int dim>
dealii::Point<dim>
SimplifiedGrainRepresentation<dim>::getCenter() const
{
  return center;
}

template <int dim>
double
SimplifiedGrainRepresentation<dim>::getRadius() const
{
  return radius;
}

template <int dim>
unsigned int
SimplifiedGrainRepresentation<dim>::getGrainId() const
{
  return grain_id;
}

template <int dim>
void
SimplifiedGrainRepresentation<dim>::setGrainId(unsigned int id)
{
  grain_id = id;
}

template <int dim>
unsigned int
SimplifiedGrainRepresentation<dim>::getOrderParameterId() const
{
  return order_parameter_id;
}

template <int dim>
void
SimplifiedGrainRepresentation<dim>::setOrderParameterId(unsigned int id)
{
  order_parameter_id = id;
}

template <int dim>
unsigned int
SimplifiedGrainRepresentation<dim>::getOldOrderParameterId() const
{
  return old_order_parameter_id;
}

template <int dim>
void
SimplifiedGrainRepresentation<dim>::setDistanceToNeighbor(double dist)
{
  distance_to_neighbor_sharing_op = dist;
}

template <int dim>
double
SimplifiedGrainRepresentation<dim>::getDistanceToNeighbor() const
{
  return distance_to_neighbor_sharing_op;
}

// ============================================================================
// Methods for SimplifiedGrainManipulator
// ============================================================================

template <int dim>
void
SimplifiedGrainManipulator<dim>::reassignGrains(
  std::vector<SimplifiedGrainRepresentation<dim>> &grain_representations,
  double                                           buffer_distance,
  std::vector<unsigned int>                       &order_parameter_id_list)
{
  for (int cycle = order_parameter_id_list.size(); cycle >= 0; cycle--)
    {
      for (unsigned int g_base = 0; g_base < grain_representations.size(); g_base++)
        {
          unsigned int order_parameter_base =
            grain_representations.at(g_base).getOrderParameterId();

          for (unsigned int g_other = 0; g_other < grain_representations.size();
               g_other++)
            {
              if (g_other != g_base)
                {
                  unsigned int order_parameter_other =
                    grain_representations.at(g_other).getOrderParameterId();

                  // Check for overlap between the base grain and the other
                  // grain
                  double center_distance =
                    grain_representations.at(g_base).getCenter().distance(
                      grain_representations.at(g_other).getCenter());
                  double sum_radii = grain_representations.at(g_base).getRadius() +
                                     grain_representations.at(g_other).getRadius();

                  if ((sum_radii + 2.0 * buffer_distance > center_distance) and
                      (order_parameter_other == order_parameter_base))
                    {
                      if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
                        {
                          std::cout << "Found overlap between grain "
                                    << grain_representations.at(g_base).getGrainId()
                                    << " and grain "
                                    << grain_representations.at(g_other).getGrainId()
                                    << " with order parameter " << order_parameter_base
                                    << std::endl;
                        }

                      grain_representations.at(g_base).setDistanceToNeighbor(
                        center_distance - sum_radii);

                      // Another loop over all of the grains to find the order
                      // parameter with the largest minimum distance to the base
                      // grain
                      std::vector<double> minimum_distance_list(
                        order_parameter_id_list.size(),
                        std::numeric_limits<double>::max());

                      for (unsigned int g_spacing_list = 0;
                           g_spacing_list < grain_representations.size();
                           g_spacing_list++)
                        {
                          if (g_spacing_list != g_base)
                            {
                              unsigned int order_parameter_spacing_list =
                                grain_representations.at(g_spacing_list)
                                  .getOrderParameterId();

                              double spacing =
                                grain_representations.at(g_base).getCenter().distance(
                                  grain_representations.at(g_spacing_list).getCenter()) -
                                grain_representations.at(g_base).getRadius() -
                                grain_representations.at(g_spacing_list).getRadius();

                              if (spacing <
                                  minimum_distance_list.at(order_parameter_spacing_list))
                                {
                                  minimum_distance_list.at(order_parameter_spacing_list) =
                                    spacing;
                                }
                            }
                        }
                      // Pick the max value of minimum_distance_list to
                      // determine which order parameter to switch the base
                      // grain to Reassign the order parameter for the grains
                      // with the conflicts with the most other order
                      // parameters. In the very last cycle, the grains that
                      // only have conflicts in their own order parameter are
                      // reassigned.
                      double       max_distance    = -std::numeric_limits<double>::max();
                      unsigned int new_op_index    = 0;
                      int          overlap_counter = 0;
                      for (unsigned int op = 0; op < minimum_distance_list.size(); op++)
                        {
                          if (minimum_distance_list.at(op) > max_distance)
                            {
                              max_distance = minimum_distance_list.at(op);
                              new_op_index = op;
                            }
                          if (minimum_distance_list.at(op) < 0)
                            {
                              overlap_counter++;
                            }
                        }
                      if (overlap_counter >= cycle)
                        {
                          grain_representations.at(g_base).setOrderParameterId(
                            new_op_index);
                          order_parameter_base = new_op_index;

                          if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                              0)
                            {
                              std::cout
                                << "Reassigning grain "
                                << grain_representations.at(g_base).getGrainId()
                                << " from order parameter "
                                << grain_representations.at(g_base)
                                     .getOldOrderParameterId()
                                << " to order parameter "
                                << grain_representations.at(g_base).getOrderParameterId()
                                << std::endl
                                << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
}

template <int dim>
void
SimplifiedGrainManipulator<dim>::transferGrainIds(
  const std::vector<SimplifiedGrainRepresentation<dim>> &old_grain_representations,
  std::vector<SimplifiedGrainRepresentation<dim>>       &new_grain_representations) const
{
  for (unsigned int g_new = 0; g_new < new_grain_representations.size(); g_new++)
    {
      double       min_distance          = std::numeric_limits<double>::max();
      unsigned int index_at_min_distance = 0;

      for (unsigned int g_old = 0; g_old < old_grain_representations.size(); g_old++)
        {
          double distance = new_grain_representations.at(g_new).getCenter().distance(
            old_grain_representations.at(g_old).getCenter());

          if (distance < min_distance)
            {
              min_distance          = distance;
              index_at_min_distance = old_grain_representations.at(g_old).getGrainId();
            }
        }
      new_grain_representations.at(g_new).setGrainId(index_at_min_distance);
    }
}

// Template instantiations
template class SimplifiedGrainManipulator<2>;
template class SimplifiedGrainManipulator<3>;

template class SimplifiedGrainRepresentation<2>;
template class SimplifiedGrainRepresentation<3>;
