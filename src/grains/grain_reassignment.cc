// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/grains/grain_reassignment.h>
#include <prismspf/grains/grains.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
void
GrainReassignmentManager<dim, degree, number>::reassign_grains(
  const SolveContext<dim, degree, number>         &solve_context,
  std::vector<SimplifiedGrainRepresentation<dim>> &simplified_grain_representations)
{
  Timer::Scope reassign_grains_scope("Grain Reassignment");

  Logger::instance() << "Reassigning grains...\n";

  // Set up references from the solve context
  const UserInputParameters<dim> &user_inputs = solve_context.get_user_inputs();

  // Create the simplified grain representations
  FloodFiller<dim, degree, number> flood_filler;
  std::vector<GrainSet<dim>>       grain_sets;

  // Search over all fields and calculate a list of GrainSet for each field assigned to
  // any grain reassignment block. One GrainSet represents one grain.
  unsigned int op_list_index             = 0;
  unsigned int last_remapped_field_index = 0;
  int          number_of_remapped_fields = 0;
  for (unsigned int field_index = 0;
       field_index < solve_context.get_field_attributes().size();
       field_index++)
    {
      const auto &variable = solve_context.get_field_attributes()[field_index];

      if (variable.grain_reassignment_block_id != -1)
        {
          number_of_remapped_fields++;
          last_remapped_field_index = field_index;

          std::vector<GrainSet<dim>> single_op_grain_sets;
          flood_filler.calc_grain_sets(
            solve_context.get_dof_manager().get_field_dof_handler(field_index),
            solve_context.get_solution_indexer().get_solution_vector(field_index),
            user_inputs.grain_reassignment_parameters.order_parameter_threshold,
            1.0 + user_inputs.grain_reassignment_parameters.order_parameter_threshold,
            0,
            field_index,
            single_op_grain_sets);

          grain_sets.insert(grain_sets.end(),
                            single_op_grain_sets.begin(),
                            single_op_grain_sets.end());
        }
    }

  DEBUG_ASSERT(number_of_remapped_fields > 0,
               "Grain reassignment active with no fields set to be remapped!!");

  // Set the grain indices to unique values
  for (unsigned int g = 0; g < grain_sets.size(); g++)
    {
      grain_sets.at(g).set_grain_index(g);
    }

  std::vector<SimplifiedGrainRepresentation<dim>> old_grain_representations =
    simplified_grain_representations;
  simplified_grain_representations.clear();
  for (unsigned int g = 0; g < grain_sets.size(); g++)
    {
      SimplifiedGrainRepresentation<dim> simplified_grain_representation(
        grain_sets.at(g));

      Logger::instance() << "Grain: " << simplified_grain_representation.get_grain_id()
                         << " "
                         << simplified_grain_representation.get_order_parameter_id()
                         << " Center: " << simplified_grain_representation.get_center()(0)
                         << " " << simplified_grain_representation.get_center()(1)
                         << std::endl;

      simplified_grain_representations.push_back(simplified_grain_representation);
    }

  if (solve_context.get_simulation_timer().get_increment() > 0 ||
      user_inputs.grain_reassignment_parameters.load_grain_structure)
    {
      SimplifiedGrainManipulator<dim>::transfer_grain_ids(
        old_grain_representations,
        simplified_grain_representations);
    }

  SimplifiedGrainManipulator<dim>::reassign_grains(
    simplified_grain_representations,
    user_inputs.grain_reassignment_parameters.exclusion_distance,
    number_of_remapped_fields,
    solve_context.get_field_attributes());

  Logger::instance() << "Grain list after reassignment:" << std::endl;

  for (unsigned int g = 0; g < simplified_grain_representations.size(); g++)
    {
      Logger::instance() << "Grain: "
                         << simplified_grain_representations[g].get_grain_id() << " "
                         << simplified_grain_representations[g].get_order_parameter_id()
                         << " Center: "
                         << simplified_grain_representations[g].get_center()(0) << " "
                         << simplified_grain_representations[g].get_center()(1)
                         << std::endl;
    }

  OrderParameterRemapper<dim, number>::remap(
    simplified_grain_representations,
    solve_context.get_solution_indexer(),
    solve_context.get_dof_manager().get_field_dof_handler(last_remapped_field_index),
    SystemWide<dim, degree>::fe_systems[0].dofs_per_cell);

  Logger::instance() << "Reassigning grains completed.\n\n";
}

template <unsigned int dim, typename number>
void
OrderParameterRemapper<dim, number>::remap(
  std::vector<SimplifiedGrainRepresentation<dim>> &grain_representations,
  SolutionIndexer<dim, number>                    &solution_indexer,
  const dealii::DoFHandler<dim>                   &dof_handler,
  unsigned int                                     dofs_per_cell)

{
  for (unsigned int g = 0; g < grain_representations.size(); g++)
    {
      if (grain_representations.at(g).get_order_parameter_id() !=
          grain_representations.at(g).get_old_order_parameter_id())
        {
          double transfer_buffer =
            std::max(0.0, grain_representations.at(g).get_distance_to_neighbor() / 2.0);

          // For now I have two loops, one where I copy the values from the old
          // order parameter to the new one and a second where I zero out the
          // old order parameter. This separation prevents writing zero-out
          // values to the new order parameter. There probably is a more
          // efficient way of doing this.
          for (const auto &dof : dof_handler.active_cell_iterators())
            {
              if (dof->is_locally_owned())
                {
                  unsigned int op_new =
                    grain_representations.at(g).get_order_parameter_id();
                  unsigned int op_old =
                    grain_representations.at(g).get_old_order_parameter_id();

                  // Check if the cell is within the simplified grain
                  // representation
                  bool in_grain = true;
                  for (unsigned int v = 0;
                       v < dealii::GeometryInfo<dim>::vertices_per_cell;
                       v++)
                    {
                      if (dof->vertex(v).distance(
                            grain_representations.at(g).get_center()) >
                          grain_representations.at(g).get_radius() + transfer_buffer)
                        {
                          in_grain = false;
                          break;
                        }
                    }

                  // If it is, move the values from the old order parameter to
                  // the new order parameter
                  if (in_grain)
                    {
                      std::vector<dealii::types::global_dof_index> dof_indices(
                        dofs_per_cell,
                        0);
                      dof->get_dof_indices(dof_indices);
                      for (const auto &index : dof_indices)
                        {
                          (solution_indexer.get_solution_vector(op_new))[index] =
                            (solution_indexer.get_solution_vector(op_old))[index];
                        }
                    }
                }
            }

          for (const auto &dof : dof_handler.active_cell_iterators())
            {
              if (dof->is_locally_owned())
                {
                  unsigned int op_old =
                    grain_representations.at(g).get_old_order_parameter_id();

                  // Check if the cell is within the simplified grain
                  // representation
                  bool in_grain = true;
                  for (unsigned int v = 0;
                       v < dealii::GeometryInfo<dim>::vertices_per_cell;
                       v++)
                    {
                      if (dof->vertex(v).distance(
                            grain_representations.at(g).get_center()) >
                          grain_representations.at(g).get_radius() + transfer_buffer)
                        {
                          in_grain = false;
                          break;
                        }
                    }

                  // If it is, set the old order parameter to zero
                  if (in_grain)
                    {
                      std::vector<dealii::types::global_dof_index> dof_indices(
                        dofs_per_cell,
                        0);
                      dof->get_dof_indices(dof_indices);

                      for (const auto &index : dof_indices)
                        {
                          (solution_indexer.get_solution_vector(op_old))[index] = 0.0;
                        }
                    }
                }
            }
        }
    }
}

// ============================================================================
// Methods for SimplifiedGrainManipulator
// ============================================================================

template <unsigned int dim>
void
SimplifiedGrainManipulator<dim>::reassign_grains(
  std::vector<SimplifiedGrainRepresentation<dim>> &grain_representations,
  double                                           buffer_distance,
  int                                              number_of_remapped_fields,
  const std::vector<FieldAttributes>              &field_attributes)
{
  for (int cycle = number_of_remapped_fields; cycle >= 0; cycle--)
    {
      for (unsigned int g_base = 0; g_base < grain_representations.size(); g_base++)
        {
          unsigned int order_parameter_base =
            grain_representations.at(g_base).get_order_parameter_id();

          for (unsigned int g_other = 0; g_other < grain_representations.size();
               g_other++)
            {
              if (g_other != g_base)
                {
                  unsigned int order_parameter_other =
                    grain_representations.at(g_other).get_order_parameter_id();

                  // Check for overlap between the base grain and the other
                  // grain
                  double center_distance =
                    grain_representations.at(g_base).get_center().distance(
                      grain_representations.at(g_other).get_center());
                  double sum_radii = grain_representations.at(g_base).get_radius() +
                                     grain_representations.at(g_other).get_radius();

                  if ((sum_radii + 2.0 * buffer_distance > center_distance) and
                      (order_parameter_other == order_parameter_base))
                    {
                      Logger::instance()
                        << "Found overlap between grain "
                        << grain_representations.at(g_base).get_grain_id()
                        << " and grain "
                        << grain_representations.at(g_other).get_grain_id()
                        << " with order parameter " << order_parameter_base << std::endl;

                      grain_representations.at(g_base).set_distance_to_neighbor(
                        center_distance - sum_radii);

                      // Another loop over all of the grains to find the order
                      // parameter with the largest minimum distance to the base
                      // grain
                      std::vector<double> minimum_distance_list(
                        number_of_remapped_fields,
                        std::numeric_limits<double>::max());

                      for (unsigned int g_spacing_list = 0;
                           g_spacing_list < grain_representations.size();
                           g_spacing_list++)
                        {
                          if (g_spacing_list != g_base)
                            {
                              unsigned int order_parameter_spacing_list =
                                grain_representations.at(g_spacing_list)
                                  .get_order_parameter_id();

                              // If this grain has a different grain_reassignment_block_id
                              // from the base grain, set the minimum distance to
                              // -double::max() such that the base grain cannot be
                              // reassigned to this OP
                              if (field_attributes.at(order_parameter_spacing_list)
                                    .grain_reassignment_block_id !=
                                  field_attributes.at(order_parameter_base)
                                    .grain_reassignment_block_id)
                                {
                                  minimum_distance_list.at(order_parameter_spacing_list) =
                                    -std::numeric_limits<double>::max();
                                }
                              else
                                {
                                  double spacing =
                                    grain_representations.at(g_base)
                                      .get_center()
                                      .distance(grain_representations.at(g_spacing_list)
                                                  .get_center()) -
                                    grain_representations.at(g_base).get_radius() -
                                    grain_representations.at(g_spacing_list).get_radius();

                                  if (spacing < minimum_distance_list.at(
                                                  order_parameter_spacing_list))
                                    {
                                      minimum_distance_list.at(
                                        order_parameter_spacing_list) = spacing;
                                    }
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
                          grain_representations.at(g_base).set_order_parameter_id(
                            new_op_index);
                          order_parameter_base = new_op_index;

                          Logger::instance()
                            << "Reassigning grain "
                            << grain_representations.at(g_base).get_grain_id()
                            << " from order parameter "
                            << grain_representations.at(g_base)
                                 .get_old_order_parameter_id()
                            << " to order parameter "
                            << grain_representations.at(g_base).get_order_parameter_id()
                            << std::endl
                            << std::endl;
                        }
                    }
                }
            }
        }
    }
}

template <unsigned int dim>
void
SimplifiedGrainManipulator<dim>::transfer_grain_ids(
  const std::vector<SimplifiedGrainRepresentation<dim>> &old_grain_representations,
  std::vector<SimplifiedGrainRepresentation<dim>>       &new_grain_representations)
{
  for (unsigned int g_new = 0; g_new < new_grain_representations.size(); g_new++)
    {
      double       min_distance          = std::numeric_limits<double>::max();
      unsigned int index_at_min_distance = 0;

      for (unsigned int g_old = 0; g_old < old_grain_representations.size(); g_old++)
        {
          double distance = new_grain_representations.at(g_new).get_center().distance(
            old_grain_representations.at(g_old).get_center());

          if (distance < min_distance)
            {
              min_distance          = distance;
              index_at_min_distance = old_grain_representations.at(g_old).get_grain_id();
            }
        }
      new_grain_representations.at(g_new).set_grain_id(index_at_min_distance);
    }
}

#include "grains/grain_reassignment.inst"

PRISMS_PF_END_NAMESPACE
