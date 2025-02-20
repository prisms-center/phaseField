// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/point.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/numerics/vector_tools_boundary.h>

#include <prismspf/config.h>
#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/nonuniform_dirichlet.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <cmath>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
constraintHandler<dim>::constraintHandler(const userInputParameters<dim> &_user_inputs)
  : user_inputs(_user_inputs)
  , constraints(user_inputs.var_attributes.size(), dealii::AffineConstraints<double>())
{}

template <int dim>
std::vector<const dealii::AffineConstraints<double> *>
constraintHandler<dim>::get_constraints()
{
  std::vector<const dealii::AffineConstraints<double> *> temp;
  temp.reserve(constraints.size());
  for (const auto &constraint : constraints)
    {
      temp.push_back(&constraint);
    }
  return temp;
}

template <int dim>
const dealii::AffineConstraints<double> &
constraintHandler<dim>::get_constraint(const unsigned int &index) const
{
  Assert(constraints.size() > index,
         dealii::ExcMessage("The constraint set does not contain index = " +
                            std::to_string(index)));
  return constraints.at(index);
}

template <int dim>
void
constraintHandler<dim>::make_constraints(
  const dealii::Mapping<dim>                   &mapping,
  const std::vector<dealii::DoFHandler<dim> *> &dof_handlers)
{
  for (const auto &[index, variable] : user_inputs.var_attributes)
    {
      make_constraint(mapping, *dof_handlers.at(index), index);
    }
}

template <int dim>
void
constraintHandler<dim>::make_constraint(const dealii::Mapping<dim>    &mapping,
                                        const dealii::DoFHandler<dim> &dof_handler,
                                        const unsigned int            &index)
{
  // Clear constraints
  constraints.at(index).clear();

  // Reinitialize constraints
  constraints.at(index).reinit(dof_handler.locally_owned_dofs(),
                               dealii::DoFTools::extract_locally_relevant_dofs(
                                 dof_handler));

  // Make hanging node constraints
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints.at(index));

  // First check the normal boundary conditions
  const auto &boundary_condition =
    user_inputs.boundary_parameters.boundary_condition_list.at(index);

  for (const auto &[component, condition] : boundary_condition)
    {
      for (const auto &[boundary_id, boundary_type] : condition.boundary_condition_map)
        {
          // Create a mask. This is only applied for vector fields to apply boundary
          // conditions to
          // each component of the vector.
          std::vector<bool> mask(dim, false);
          mask.at(component) = true;

          if (boundary_type == boundaryCondition::type::NATURAL)
            {
              // Do nothing because they are naturally enforced.
              continue;
            }
          else if (boundary_type == boundaryCondition::type::DIRICHLET)
            {
              if (user_inputs.var_attributes.at(index).field_type != fieldType::VECTOR)
                {
                  dealii::VectorTools::interpolate_boundary_values(
                    mapping,
                    dof_handler,
                    boundary_id,
                    dealii::Functions::ConstantFunction<dim>(
                      condition.dirichlet_value_map.at(boundary_id),
                      1),
                    constraints.at(index));
                }
              else
                {
                  dealii::VectorTools::interpolate_boundary_values(
                    mapping,
                    dof_handler,
                    boundary_id,
                    dealii::Functions::ConstantFunction<dim>(
                      condition.dirichlet_value_map.at(boundary_id),
                      dim),
                    constraints.at(index),
                    mask);
                }
            }
          else if (boundary_type == boundaryCondition::type::PERIODIC)
            {
              // Skip boundary ids that are odd since those map to the even faces
              if (boundary_id % 2 != 0)
                {
                  continue;
                }
              // Create a vector of matched pairs that we fill and enforce upon the
              // constaints
              std::vector<dealii::GridTools::PeriodicFacePair<
                typename dealii::DoFHandler<dim>::cell_iterator>>
                periodicity_vector;

              // Determine the direction
              const auto direction =
                static_cast<unsigned int>(std::floor(boundary_id / dim));

              // Collect the matched pairs on the coarsest level of the mesh
              dealii::GridTools::collect_periodic_faces(dof_handler,
                                                        boundary_id,
                                                        boundary_id + 1,
                                                        direction,
                                                        periodicity_vector);

              // Set constraints
              if (user_inputs.var_attributes.at(index).field_type != fieldType::VECTOR)
                {
                  dealii::DoFTools::make_periodicity_constraints<dim, dim>(
                    periodicity_vector,
                    constraints.at(index));
                }
              else
                {
                  dealii::DoFTools::make_periodicity_constraints<dim, dim>(
                    periodicity_vector,
                    constraints.at(index),
                    mask);
                }
            }
          else if (boundary_type == boundaryCondition::type::NEUMANN)
            {
              Assert(false, FeatureNotImplemented("Neumann boundary conditions"));
            }
          else if (boundary_type == boundaryCondition::type::NON_UNIFORM_DIRICHLET)
            {
              if (user_inputs.var_attributes.at(index).field_type != fieldType::VECTOR)
                {
                  dealii::VectorTools::interpolate_boundary_values(
                    mapping,
                    dof_handler,
                    boundary_id,
                    nonuniformDirichlet<dim, fieldType::SCALAR>(index, boundary_id),
                    constraints.at(index));
                }
              else
                {
                  dealii::VectorTools::interpolate_boundary_values(
                    mapping,
                    dof_handler,
                    boundary_id,
                    nonuniformDirichlet<dim, fieldType::VECTOR>(index, boundary_id),
                    constraints.at(index),
                    mask);
                }
            }
          else if (boundary_type == boundaryCondition::type::NON_UNIFORM_NEUMANN)
            {
              Assert(false,
                     FeatureNotImplemented("Nonuniform neumann boundary conditions"));
            }
        }
    }

  // Second check for pinned points, if they exist
  if (user_inputs.boundary_parameters.pinned_point_list.find(index) !=
      user_inputs.boundary_parameters.pinned_point_list.end())
    {
      set_pinned_point(dof_handler, index);
    }

  // Close constraints
  constraints.at(index).close();
}

template <int dim>
void
constraintHandler<dim>::set_pinned_point(const dealii::DoFHandler<dim> &dof_handler,
                                         const unsigned int            &index)
{
  Assert(user_inputs.var_attributes.at(index).field_type == fieldType::VECTOR,
         FeatureNotImplemented("Pinned points for vector fields"));

  const auto &value_point_pair =
    user_inputs.boundary_parameters.pinned_point_list.at(index);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int i = 0; i < dealii::GeometryInfo<dim>::vertices_per_cell; ++i)
            {
              // Check if the vertex is the target vertex
              if (value_point_pair.second.distance(cell->vertex(i)) <
                  1.0e-2 * cell->diameter())
                {
                  unsigned int nodeID = cell->vertex_dof_index(i, 0);
                  constraints.at(index).add_line(nodeID);
                  constraints.at(index).set_inhomogeneity(nodeID, value_point_pair.first);
                }
            }
        }
    }
}

INSTANTIATE_UNI_TEMPLATE(constraintHandler)

PRISMS_PF_END_NAMESPACE