#include "../../include/OrderParameterRemapper.h"

template <int dim>
void
OrderParameterRemapper<dim>::remap(
  std::vector<SimplifiedGrainRepresentation<dim>> &grain_representations,
  std::vector<vectorType *>                       &solution_fields,
  dealii::DoFHandler<dim>                         &dof_handler,
  unsigned int                                     dofs_per_cell,
  double                                           buffer)
{
  for (unsigned int g = 0; g < grain_representations.size(); g++)
    {
      if (grain_representations.at(g).getOrderParameterId() !=
          grain_representations.at(g).getOldOrderParameterId())
        {
          double transfer_buffer =
            std::max(0.0, grain_representations.at(g).getDistanceToNeighbor() / 2.0);

          typename dealii::DoFHandler<dim>::active_cell_iterator di =
            dof_handler.begin_active();

          // For now I have two loops, one where I copy the values from the old
          // order parameter to the new one and a second where I zero out the
          // old order parameter. This separation prevents writing zero-out
          // values to the new order parameter. There probably is a more
          // efficient way of doing this.
          while (di != dof_handler.end())
            {
              if (di->is_locally_owned())
                {
                  unsigned int op_new = grain_representations.at(g).getOrderParameterId();
                  unsigned int op_old =
                    grain_representations.at(g).getOldOrderParameterId();

                  // Check if the cell is within the simplified grain
                  // representation
                  bool in_grain = true;
                  for (unsigned int v = 0;
                       v < dealii::GeometryInfo<dim>::vertices_per_cell;
                       v++)
                    {
                      if (di->vertex(v).distance(
                            grain_representations.at(g).getCenter()) >
                          grain_representations.at(g).getRadius() + transfer_buffer)
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
                      di->get_dof_indices(dof_indices);
                      for (unsigned int i = 0; i < dof_indices.size(); i++)
                        {
                          (*solution_fields.at(op_new))[dof_indices.at(i)] =
                            (*solution_fields.at(op_old))[dof_indices.at(i)];
                        }
                    }
                }
              ++di;
            }

          di = dof_handler.begin_active();

          while (di != dof_handler.end())
            {
              if (di->is_locally_owned())
                {
                  unsigned int op_new = grain_representations.at(g).getOrderParameterId();
                  unsigned int op_old =
                    grain_representations.at(g).getOldOrderParameterId();

                  // Check if the cell is within the simplified grain
                  // representation
                  bool in_grain = true;
                  for (unsigned int v = 0;
                       v < dealii::GeometryInfo<dim>::vertices_per_cell;
                       v++)
                    {
                      if (di->vertex(v).distance(
                            grain_representations.at(g).getCenter()) >
                          grain_representations.at(g).getRadius() + transfer_buffer)
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
                      di->get_dof_indices(dof_indices);

                      for (unsigned int i = 0; i < dof_indices.size(); i++)
                        {
                          (*solution_fields.at(op_old))[dof_indices.at(i)] = 0.0;
                        }
                    }
                }
              ++di;
            }
        }
    }
}

template <int dim>
void
OrderParameterRemapper<dim>::remap_from_index_field(
  std::vector<SimplifiedGrainRepresentation<dim>> &grain_representations,
  const vectorType                                *grain_index_field,
  std::vector<vectorType *>                       &solution_fields,
  dealii::DoFHandler<dim>                         &dof_handler,
  unsigned int                                     dofs_per_cell,
  double                                           buffer)
{
  for (unsigned int g = 0; g < grain_representations.size(); g++)
    {
      std::cout << "Grain: " << grain_representations.at(g).getGrainId()
                << " Old OP: " << grain_representations.at(g).getOldOrderParameterId()
                << " New OP: " << grain_representations.at(g).getOrderParameterId()
                << std::endl;

      double transfer_buffer =
        std::max(0.0, grain_representations.at(g).getDistanceToNeighbor() / 2.0);

      typename dealii::DoFHandler<dim>::active_cell_iterator di =
        dof_handler.begin_active();

      // For now I have two loops, one where I copy the values from the old
      // order parameter to the new one and a second where I zero out the old
      // order parameter. This separation prevents writing zero-out values to
      // the new order parameter. There probably is a more efficient way of
      // doing this.
      while (di != dof_handler.end())
        {
          if (di->is_locally_owned())
            {
              unsigned int op_new = grain_representations.at(g).getOrderParameterId();
              unsigned int op_old = grain_representations.at(g).getOldOrderParameterId();

              // Check if the cell is within the simplified grain representation
              bool in_grain = true;
              for (unsigned int v = 0; v < dealii::GeometryInfo<dim>::vertices_per_cell;
                   v++)
                {
                  if (di->vertex(v).distance(grain_representations.at(g).getCenter()) >
                      grain_representations.at(g).getRadius() + transfer_buffer)
                    {
                      in_grain = false;
                      break;
                    }
                }

              // If it is, move the values from the old order parameter to the
              // new order parameter
              if (in_grain)
                {
                  std::vector<dealii::types::global_dof_index> dof_indices(dofs_per_cell,
                                                                           0);
                  di->get_dof_indices(dof_indices);
                  for (unsigned int i = 0; i < dof_indices.size(); i++)
                    {
                      if (std::abs((*grain_index_field)[dof_indices.at(i)] -
                                   (double) grain_representations.at(g).getGrainId()) <
                          1e-6)
                        {
                          (*solution_fields.at(op_new))[dof_indices.at(i)] = 1.0;
                        }
                    }
                }
            }
          ++di;
        }
    }
}

template class OrderParameterRemapper<2>;
template class OrderParameterRemapper<3>;
