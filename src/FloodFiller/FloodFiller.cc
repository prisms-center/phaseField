#include "../../include/FloodFiller.h"

#include <numeric>

template <int dim, int degree>
void
FloodFiller<dim, degree>::calcGrainSets(
  [[maybe_unused]] dealii::FESystem<dim> &finite_element,
  dealii::DoFHandler<dim>                &dof_handler,
  vectorType                             *solution_field,
  double                                  threshold_lower,
  double                                  threshold_upper,
  int                                     min_id,
  unsigned int                            order_parameter_index,
  std::vector<GrainSet<dim>>             &grain_sets)
{
  unsigned int grain_index = 0;

  // Loop through the whole mesh and set the user flags to false (so everything
  // is considered unmarked)
  typename dealii::DoFHandler<dim>::cell_iterator cell = dof_handler.begin();
  while (cell != dof_handler.end())
    {
      cell->clear_user_flag();
      ++cell;
    }

  GrainSet<dim> grain_set;
  grain_sets.push_back(grain_set);
  grain_sets.back().setOrderParameterIndex(order_parameter_index);

  // The flood fill loop
  cell = dof_handler.begin();
  while (cell != dof_handler.end())
    {
      if (!cell->has_children())
        {
          bool grain_assigned = false;
          recursiveFloodFill<typename dealii::DoFHandler<dim>::cell_iterator>(
            cell,
            dof_handler.end(),
            solution_field,
            threshold_lower,
            threshold_upper,
            min_id,
            grain_index,
            grain_sets,
            grain_assigned);

          if (grain_assigned)
            {
              // Get the grain set initialized for the next grain to be found
              grain_index++;
              GrainSet<dim> new_grain_set;
              new_grain_set.setOrderParameterIndex(order_parameter_index);
              grain_sets.push_back(new_grain_set);
            }
        }

      ++cell;
    }

  // If the last grain was initialized but empty, delete it
  if (grain_sets.back().getVertexList().size() == 0)
    {
      grain_sets.pop_back();
    }

  // Generate global list of the grains & send the grain set info to all processors so
  // everyone has the full list
  if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) > 1)
    {
      createGlobalGrainSetList(grain_sets);
    }

  // Merge grains sharing common vertices
  mergeSplitGrains(grain_sets);
}

// NOLINTBEGIN(misc-no-recursion)
template <int dim, int degree>
template <typename T>
void
FloodFiller<dim, degree>::recursiveFloodFill(T                           cell,
                                             T                           cell_end,
                                             vectorType                 *solution_field,
                                             double                      threshold_lower,
                                             double                      threshold_upper,
                                             int                         min_id,
                                             unsigned int               &grain_index,
                                             std::vector<GrainSet<dim>> &grain_sets,
                                             bool                       &grain_assigned)
{
  if (cell == cell_end)
    {
      return;
    }

  // Check if the cell has been marked yet
  if (cell->user_flag_set())
    {
      return;
    }

  if (cell->has_children())
    {
      // Call recursiveFloodFill on the element's children
      for (unsigned int n_child = 0; n_child < cell->n_children(); n_child++)
        {
          recursiveFloodFill<T>(cell->child(n_child),
                                cell_end,
                                solution_field,
                                threshold_lower,
                                threshold_upper,
                                min_id,
                                grain_index,
                                grain_sets,
                                grain_assigned);
        }
    }
  else
    {
      if (!cell->is_locally_owned())
        {
          return;
        }

      cell->set_user_flag();

      dealii::FEValues<dim> fe_values(*fe, quadrature, dealii::update_values);

      std::vector<double>             var_values(num_quad_points);
      std::vector<dealii::Point<dim>> q_point_list(num_quad_points);

      // Get the most common value for the element
      fe_values.reinit(cell);
      fe_values.get_function_values(*solution_field, var_values);

      std::map<double, int> quadratureValues;
      int                   maxNumberSeen         = 0;
      double                mostCommonQPointValue = -1;
      for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
        {
          // Add the number of times that var_values[q_point] has
          // been seen
          if (var_values[q_point] > min_id)
            {
              ++quadratureValues[var_values[q_point]];
            }
          if (quadratureValues[var_values[q_point]] > maxNumberSeen)
            {
              maxNumberSeen         = quadratureValues[var_values[q_point]];
              mostCommonQPointValue = var_values[q_point];
            }
        }
      double ele_val = mostCommonQPointValue;

      if (ele_val > threshold_lower && ele_val < threshold_upper)
        {
          grain_assigned = true;

          std::vector<dealii::Point<dim>> vertex_list;
          for (unsigned int vertex_index = 0;
               vertex_index < dealii::Utilities::fixed_power<dim>(2.0);
               vertex_index++)
            {
              vertex_list.push_back(cell->vertex(vertex_index));
            }
          grain_sets.back().addVertexList(vertex_list);

          // Call recursiveFloodFill on the element's neighbors
          for (unsigned int n_child = 0; n_child < 2 * dim; n_child++)
            {
              recursiveFloodFill<T>(cell->neighbor(n_child),
                                    cell_end,
                                    solution_field,
                                    threshold_lower,
                                    threshold_upper,
                                    min_id,
                                    grain_index,
                                    grain_sets,
                                    grain_assigned);
            }
        }
    }
}

// NOLINTEND(misc-no-recursion)

// =================================================================================
// All-to-all communication of the grain sets
// =================================================================================
template <int dim, int degree>
void
FloodFiller<dim, degree>::createGlobalGrainSetList(
  std::vector<GrainSet<dim>> &grain_sets) const
{
  unsigned int numProcs = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  unsigned int num_grains_local = grain_sets.size();

  // Convert the grain_set object into a group of vectors
  std::vector<unsigned int> order_parameters;
  std::vector<unsigned int> num_elements;
  std::vector<double>       vertices;

  for (unsigned int g = 0; g < grain_sets.size(); g++)
    {
      order_parameters.push_back(grain_sets.at(g).getOrderParameterIndex());

      std::vector<std::vector<dealii::Point<dim>>> vertex_list =
        grain_sets[g].getVertexList();
      num_elements.push_back(vertex_list.size());

      for (unsigned int c = 0; c < num_elements[g]; c++)
        {
          for (unsigned int v = 0; v < dealii::Utilities::fixed_power<dim>(2.0); v++)
            {
              for (unsigned int d = 0; d < dim; d++)
                {
                  vertices.push_back(vertex_list[c][v][d]);
                }
            }
        }
    }

  unsigned int num_vertices = 0;
  for (unsigned int g = 0; g < grain_sets.size(); g++)
    {
      num_vertices += num_elements[g] * dealii::Utilities::fixed_power<dim>(2) * dim;
    }

  // Communicate how many grains each core has
  std::vector<int> num_grains_per_core(numProcs, 0);

  MPI_Allgather(&num_grains_local,
                1,
                MPI_INT,
                num_grains_per_core.data(),
                1,
                MPI_INT,
                MPI_COMM_WORLD);

  int num_grains_global =
    std::accumulate(num_grains_per_core.begin(), num_grains_per_core.end(), 0);

  // Communicate the order_parameters
  std::vector<int> offset(numProcs, 0);
  for (int n = 1; n < numProcs; n++)
    {
      offset[n] = offset[n - 1] + num_grains_per_core[n - 1];
    }

  std::vector<unsigned int> order_parameters_global(num_grains_global, 0);

  MPI_Allgatherv(order_parameters.data(),
                 num_grains_local,
                 MPI_UNSIGNED,
                 order_parameters_global.data(),
                 num_grains_per_core.data(),
                 offset.data(),
                 MPI_UNSIGNED,
                 MPI_COMM_WORLD);

  // Communicate the number of elements
  std::vector<unsigned int> num_elements_global(num_grains_global, 0);

  MPI_Allgatherv(num_elements.data(),
                 num_grains_local,
                 MPI_UNSIGNED,
                 num_elements_global.data(),
                 num_grains_per_core.data(),
                 offset.data(),
                 MPI_UNSIGNED,
                 MPI_COMM_WORLD);

  // Communicate the vertices
  unsigned int total_elements =
    std::accumulate(num_elements_global.begin(), num_elements_global.end(), 0);
  int num_vertices_global = (unsigned int) total_elements *
                            dealii::Utilities::fixed_power<dim>(2) * (unsigned int) dim;
  std::vector<double> vertices_global(num_vertices_global, 0);

  std::vector<int> num_vertices_per_core;

  unsigned int g = 0;
  for (unsigned int i = 0; i < numProcs; i++)
    {
      int num_vert_single_core = 0;
      for (int j = 0; j < num_grains_per_core.at(i); j++)
        {
          num_vert_single_core += num_elements_global.at(g) *
                                  dealii::Utilities::fixed_power<dim>(2) *
                                  (unsigned int) dim;
          g++;
        }
      num_vertices_per_core.push_back(num_vert_single_core);
    }

  offset.at(0) = 0;
  for (unsigned int n = 1; n < numProcs; n++)
    {
      offset[n] = offset[n - 1] + num_vertices_per_core[n - 1];
    }

  MPI_Allgatherv(vertices.data(),
                 num_vertices,
                 MPI_DOUBLE,
                 vertices_global.data(),
                 num_vertices_per_core.data(),
                 offset.data(),
                 MPI_DOUBLE,
                 MPI_COMM_WORLD);

  // Put the GrainSet objects back together
  grain_sets.clear();

  for (int g = 0; g < num_grains_global; g++)
    {
      GrainSet<dim> new_grain_set;
      for (unsigned int c = 0; c < num_elements_global.at(g); c++)
        {
          std::vector<dealii::Point<dim>> verts;
          for (unsigned int v = 0; v < dealii::Utilities::fixed_power<dim>(2.0); v++)
            {
              double coords[dim];
              for (unsigned int d = 0; d < dim; d++)
                {
                  coords[d] = vertices_global.front();
                  vertices_global.erase(vertices_global.begin());
                }
              dealii::Tensor<1, dim> tensor_coords(coords);
              dealii::Point<dim>     vert(tensor_coords);
              verts.push_back(vert);
            }
          new_grain_set.addVertexList(verts);
        }
      new_grain_set.setOrderParameterIndex(order_parameters_global.at(g));
      grain_sets.push_back(new_grain_set);
    }
}

// =================================================================================
// Check to see if any grains on different processors share vertices
// =================================================================================

template <int dim, int degree>
void
FloodFiller<dim, degree>::mergeSplitGrains(std::vector<GrainSet<dim>> &grain_sets) const
{
  // Loop though each vertex in the base grain "g"
  for (unsigned int g = 0; g < grain_sets.size(); g++)
    {
      std::vector<std::vector<dealii::Point<dim>>> vertex_list =
        grain_sets[g].getVertexList();

      // Now cycle through the other grains to find overlapping elements
      for (unsigned int g_other = g + 1; g_other < grain_sets.size(); g_other++)
        {
          bool matching_vert = false;

          std::vector<std::vector<dealii::Point<dim>>> vertex_list_other =
            grain_sets[g_other].getVertexList();

          for (unsigned int c = 0; c < vertex_list.size(); c++)
            {
              for (unsigned int v = 0; v < dealii::Utilities::fixed_power<dim>(2.0); v++)
                {
                  for (unsigned int c_other = 0; c_other < vertex_list_other.size();
                       c_other++)
                    {
                      for (unsigned int v_other = 0;
                           v_other < dealii::Utilities::fixed_power<dim>(2.0);
                           v_other++)
                        {
                          // Check if the vertices match
                          if (vertex_list[c][v] == vertex_list_other[c_other][v_other])
                            {
                              matching_vert = true;
                              break;
                            }
                          if (matching_vert)
                            {
                              break;
                            }
                        }
                      if (matching_vert)
                        {
                          break;
                        }
                    }
                  if (matching_vert)
                    {
                      break;
                    }
                }
              if (matching_vert)
                {
                  break;
                }
            }

          if (matching_vert)
            {
              for (unsigned int c_base = 0; c_base < vertex_list.size(); c_base++)
                {
                  grain_sets[g_other].addVertexList(vertex_list.at(c_base));
                }
              grain_sets.erase(grain_sets.begin() + g);
              g--;
              break;
            }
        }
    }
}

// Template instantiations
template class FloodFiller<2, 1>;
template class FloodFiller<3, 1>;

template class FloodFiller<2, 2>;
template class FloodFiller<3, 2>;

template class FloodFiller<2, 3>;
template class FloodFiller<3, 3>;

template class FloodFiller<2, 4>;
template class FloodFiller<3, 4>;

template class FloodFiller<2, 5>;
template class FloodFiller<3, 5>;

template class FloodFiller<2, 6>;
template class FloodFiller<3, 6>;
