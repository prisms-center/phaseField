// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/numerics/vector_tools_boundary.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/nonuniform_dirichlet.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/boundary_parameters.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <algorithm>
#include <iterator>
#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
constraintHandler<dim>::constraintHandler(const userInputParameters<dim> &_user_inputs,
                                          const MGInfo<dim>              &_mg_info)
  : user_inputs(&_user_inputs)
  , mg_info(&_mg_info)
  , constraints(_user_inputs.var_attributes->size(), dealii::AffineConstraints<double>())
{
  // If we don't have multigrid, we can return early
  if (!_mg_info.has_multigrid())
    {
      return;
    }
  has_multigrid = true;

  global_min_level = _mg_info.get_mg_min_level();
  mg_constraints.resize(_mg_info.get_mg_depth());
  for (unsigned int level = 0; level < _mg_info.get_mg_depth(); level++)
    {
      mg_constraints[level].resize(_mg_info.get_mg_breadth(level));
    }
}

template <int dim>
std::vector<const dealii::AffineConstraints<double> *>
constraintHandler<dim>::get_constraints()
{
  std::vector<const dealii::AffineConstraints<double> *> temp;
  temp.reserve(constraints.size());

  std::transform(constraints.begin(),
                 constraints.end(),
                 std::back_inserter(temp),
                 [](const dealii::AffineConstraints<double> &constraint)
                 {
                   return &constraint;
                 });

  return temp;
}

template <int dim>
const dealii::AffineConstraints<double> &
constraintHandler<dim>::get_constraint(unsigned int index) const
{
  Assert(constraints.size() > index,
         dealii::ExcMessage("The constraint set does not contain index = " +
                            std::to_string(index)));
  return constraints.at(index);
}

template <int dim>
std::vector<const dealii::AffineConstraints<float> *>
constraintHandler<dim>::get_mg_constraints(unsigned int level)
{
  Assert(has_multigrid, dealii::ExcNotInitialized());

  // Convert the absolute level to a relative level.
  const unsigned int relative_level = level - global_min_level;
  Assert(relative_level < mg_constraints.size(),
         dealii::ExcMessage("The mg constraint set does not contain level = " +
                            std::to_string(level)));

  std::vector<const dealii::AffineConstraints<float> *> temp;
  temp.reserve(mg_constraints[relative_level].size());
  std::transform(mg_constraints[relative_level].begin(),
                 mg_constraints[relative_level].end(),
                 std::back_inserter(temp),
                 [](const dealii::AffineConstraints<float> &constraint)
                 {
                   return &constraint;
                 });

  return temp;
}

template <int dim>
const dealii::AffineConstraints<float> &
constraintHandler<dim>::get_mg_constraint(unsigned int level, unsigned int index) const
{
  Assert(has_multigrid, dealii::ExcNotInitialized());
  // Convert the absolute level to a relative level.
  const unsigned int relative_level = level - global_min_level;
  Assert(relative_level < mg_constraints.size(),
         dealii::ExcMessage("The mg constraint set does not contain level = " +
                            std::to_string(level)));
  const unsigned int global_index = mg_info->get_global_index(index, relative_level);
  Assert(index < mg_constraints[relative_level].size(),
         dealii::ExcMessage("The mg constraint set does not contain index = " +
                            std::to_string(global_index)));
  return mg_constraints[relative_level][index];
}

template <int dim>
void
constraintHandler<dim>::make_constraints(
  const dealii::Mapping<dim>                         &mapping,
  const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers)
{
  for (const auto &[index, variable] : *user_inputs->var_attributes)
    {
      make_constraint(mapping, *dof_handlers.at(index), index);
    }
}

template <int dim>
void
constraintHandler<dim>::make_mg_constraints(
  const dealii::Mapping<dim>                         &mapping,
  const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers,
  unsigned int                                        level)
{
  Assert(has_multigrid, dealii::ExcNotInitialized());

  for (unsigned int index = 0; index < dof_handlers.size(); index++)
    {
      // TODO (landinjm): Fix this so we can actually apply constraints to the LHS fields
      // baseded on whether it is normal or change. For now, I know this will have no
      // effect on the precipitate app, so I'll leave it.
      make_mg_constraint(mapping,
                         *dof_handlers[index],
                         index,
                         level,
                         dependencyType::CHANGE);
    }
}

template <int dim>
void
constraintHandler<dim>::make_constraint(const dealii::Mapping<dim>    &mapping,
                                        const dealii::DoFHandler<dim> &dof_handler,
                                        unsigned int                   index)
{
  // The constraint that we are going to modify
  Assert(index < constraints.size(),
         dealii::ExcMessage("The constraint set does not contain index = " +
                            std::to_string(index)));
  auto &local_constraint = constraints[index];

  apply_generic_constraints<double>(dof_handler, local_constraint);

  // First check the normal boundary conditions
  const auto &boundary_condition =
    user_inputs->boundary_parameters.boundary_condition_list.at(index);
  for (const auto &[component, condition] : boundary_condition)
    {
      for (const auto &[boundary_id, boundary_type] : condition.boundary_condition_map)
        {
          if (user_inputs->var_attributes->at(index).field_type != fieldType::VECTOR)
            {
              apply_constraints<double, 1>(mapping,
                                           dof_handler,
                                           local_constraint,
                                           condition,
                                           boundary_type,
                                           boundary_id,
                                           component,
                                           index);
              continue;
            }
          apply_constraints<double, dim>(mapping,
                                         dof_handler,
                                         local_constraint,
                                         condition,
                                         boundary_type,
                                         boundary_id,
                                         component,
                                         index);
        }
    }

  // Second check for pinned points, if they exist
  if (user_inputs->boundary_parameters.pinned_point_list.contains(index))
    {
      set_pinned_point<double>(dof_handler, local_constraint, index);
    }

  // Close constraints
  local_constraint.close();
}

template <int dim>
void
constraintHandler<dim>::make_mg_constraint(const dealii::Mapping<dim>    &mapping,
                                           const dealii::DoFHandler<dim> &dof_handler,
                                           unsigned int                   index,
                                           unsigned int                   level,
                                           dependencyType                 dependency_type)
{
  // Convert the global index and absolute level to a local index and relative level.
  const unsigned int relative_level = level - global_min_level;
  const unsigned int global_index   = mg_info->get_global_index(index, relative_level);

  // The constraint that we are going to modify
  Assert(relative_level < mg_constraints.size(),
         dealii::ExcMessage("The mg constraint set does not contain level = " +
                            std::to_string(level)));
  Assert(index < mg_constraints[relative_level].size(),
         dealii::ExcMessage("The mg constraint set does not contain index = " +
                            std::to_string(global_index)));
  auto &local_constraint = mg_constraints[relative_level][index];

  apply_generic_constraints<float>(dof_handler, local_constraint);

  if (dependency_type == dependencyType::CHANGE)
    {
      const auto &boundary_condition =
        this->user_inputs->boundary_parameters.boundary_condition_list.at(global_index);
      for (const auto &[component, condition] : boundary_condition)
        {
          for (const auto &[boundary_id, boundary_type] :
               condition.boundary_condition_map)
            {
              for (const auto &[boundary_id, boundary_type] :
                   condition.boundary_condition_map)
                {
                  if (user_inputs->var_attributes->at(global_index).field_type !=
                      fieldType::VECTOR)
                    {
                      apply_mg_constraints<float, 1>(mapping,
                                                     dof_handler,
                                                     local_constraint,
                                                     boundary_type,
                                                     boundary_id,
                                                     component);
                      continue;
                    }
                  apply_mg_constraints<float, dim>(mapping,
                                                   dof_handler,
                                                   local_constraint,
                                                   boundary_type,
                                                   boundary_id,
                                                   component);
                }
            }
        }

      // Second check for pinned points, if they exist
      if (user_inputs->boundary_parameters.pinned_point_list.contains(global_index))
        {
          set_mg_pinned_point<float>(dof_handler, local_constraint, global_index);
        }
    }
  else if (dependency_type == dependencyType::NORMAL)
    {
      const auto &boundary_condition =
        this->user_inputs->boundary_parameters.boundary_condition_list.at(global_index);
      for (const auto &[component, condition] : boundary_condition)
        {
          for (const auto &[boundary_id, boundary_type] :
               condition.boundary_condition_map)
            {
              for (const auto &[boundary_id, boundary_type] :
                   condition.boundary_condition_map)
                {
                  if (user_inputs->var_attributes->at(global_index).field_type !=
                      fieldType::VECTOR)
                    {
                      apply_constraints<float, 1>(mapping,
                                                  dof_handler,
                                                  local_constraint,
                                                  condition,
                                                  boundary_type,
                                                  boundary_id,
                                                  component,
                                                  index);
                      continue;
                    }
                  apply_constraints<float, dim>(mapping,
                                                dof_handler,
                                                local_constraint,
                                                condition,
                                                boundary_type,
                                                boundary_id,
                                                component,
                                                index);
                }
            }
        }

      // Second check for pinned points, if they exist
      if (user_inputs->boundary_parameters.pinned_point_list.contains(global_index))
        {
          set_pinned_point<float>(dof_handler, local_constraint, global_index);
        }
    }
  else
    {
      AssertThrow(false, dealii::ExcNotImplemented());
    }

  local_constraint.close();
}

template <int dim>
template <typename number>
void
constraintHandler<dim>::set_pinned_point(const dealii::DoFHandler<dim>     &dof_handler,
                                         dealii::AffineConstraints<number> &constraints,
                                         unsigned int                       index) const
{
  Assert(user_inputs->var_attributes->at(index).field_type == fieldType::VECTOR,
         FeatureNotImplemented("Pinned points for vector fields"));

  const number tolerance = 1.0e-2;

  const auto &value_point_pair =
    user_inputs->boundary_parameters.pinned_point_list.at(index);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int i = 0; i < dealii::GeometryInfo<dim>::vertices_per_cell; ++i)
            {
              // Check if the vertex is the target vertex
              if (value_point_pair.second.distance(cell->vertex(i)) <
                  tolerance * cell->diameter())
                {
                  const unsigned int nodeID = cell->vertex_dof_index(i, 0);
                  constraints.add_line(nodeID);
                  constraints.set_inhomogeneity(nodeID, value_point_pair.first);
                }
            }
        }
    }
}

template <int dim>
template <typename number>
void
constraintHandler<dim>::set_mg_pinned_point(
  const dealii::DoFHandler<dim>     &dof_handler,
  dealii::AffineConstraints<number> &constraints,
  unsigned int                       index) const
{
  Assert(user_inputs->var_attributes->at(index).field_type == fieldType::VECTOR,
         FeatureNotImplemented("Pinned points for vector fields"));

  const number tolerance = 1.0e-2;

  const auto &value_point_pair =
    user_inputs->boundary_parameters.pinned_point_list.at(index);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int i = 0; i < dealii::GeometryInfo<dim>::vertices_per_cell; ++i)
            {
              // Check if the vertex is the target vertex
              if (value_point_pair.second.distance(cell->vertex(i)) <
                  tolerance * cell->diameter())
                {
                  const unsigned int nodeID = cell->vertex_dof_index(i, 0);
                  constraints.add_line(nodeID);
                  constraints.set_inhomogeneity(nodeID, 0.0);
                }
            }
        }
    }
}

template <int dim>
template <typename number>
void
constraintHandler<dim>::apply_generic_constraints(
  const dealii::DoFHandler<dim>     &dof_handler,
  dealii::AffineConstraints<number> &constraints) const
{
  // Clear constraints
  constraints.clear();

  // Reinitialize constraints
  constraints.reinit(dof_handler.locally_owned_dofs(),
                     dealii::DoFTools::extract_locally_relevant_dofs(dof_handler));

  // Make hanging node constraints
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints);
}

template <int dim>
template <typename number, int spacedim>
void
constraintHandler<dim>::apply_constraints(const dealii::Mapping<dim>        &mapping,
                                          const dealii::DoFHandler<dim>     &dof_handler,
                                          dealii::AffineConstraints<number> &constraints,
                                          const boundaryCondition &boundary_condition,
                                          boundaryCondition::type  boundary_type,
                                          unsigned int             boundary_id,
                                          unsigned int             component,
                                          unsigned int             index) const
{
  constexpr bool is_vector_field = spacedim != 1;

  // Create a component mask. This will only select a certain component for vector
  // fields
  dealii::ComponentMask mask = {};
  if constexpr (is_vector_field)
    {
      std::vector<bool> temp_mask(dim, false);
      temp_mask.at(component) = true;

      mask = dealii::ComponentMask(temp_mask);
    }

  // Apply the boundary conditions
  if (boundary_type == boundaryCondition::type::NATURAL)
    {
      // Do nothing because they are naturally enforced.
      return;
    }
  if (boundary_type == boundaryCondition::type::DIRICHLET)
    {
      dealii::VectorTools::interpolate_boundary_values(
        mapping,
        dof_handler,
        boundary_id,
        dealii::Functions::ConstantFunction<dim, number>(
          boundary_condition.dirichlet_value_map.at(boundary_id),
          is_vector_field ? dim : 1),
        constraints,
        mask);
      return;
    }
  if (boundary_type == boundaryCondition::type::PERIODIC)
    {
      // Skip boundary ids that are odd since those map to the even faces
      if (boundary_id % 2 != 0)
        {
          return;
        }
      // Create a vector of matched pairs that we fill and enforce upon the
      // constaints
      std::vector<dealii::GridTools::PeriodicFacePair<
        typename dealii::DoFHandler<dim>::cell_iterator>>
        periodicity_vector;

      // Determine the direction
      const auto direction = static_cast<unsigned int>(std::floor(boundary_id / dim));

      // Collect the matched pairs on the coarsest level of the mesh
      dealii::GridTools::collect_periodic_faces(dof_handler,
                                                boundary_id,
                                                boundary_id + 1,
                                                direction,
                                                periodicity_vector);

      // Set constraints
      dealii::DoFTools::make_periodicity_constraints<dim, dim, number>(periodicity_vector,
                                                                       constraints,
                                                                       mask);
      return;
    }
  if (boundary_type == boundaryCondition::type::NEUMANN)
    {
      Assert(false, FeatureNotImplemented("Neumann boundary conditions"));
      return;
    }
  if (boundary_type == boundaryCondition::type::NON_UNIFORM_DIRICHLET)
    {
      dealii::VectorTools::interpolate_boundary_values(
        mapping,
        dof_handler,
        boundary_id,
        nonuniformDirichlet<dim, number>(index,
                                         boundary_id,
                                         *user_inputs,
                                         is_vector_field ? dim : 1),
        constraints,
        mask);
      return;
    }
  if (boundary_type == boundaryCondition::type::NON_UNIFORM_NEUMANN)
    {
      Assert(false, FeatureNotImplemented("Nonuniform neumann boundary conditions"));
      return;
    }
}

template <int dim>
template <typename number, int spacedim>
void
constraintHandler<dim>::apply_mg_constraints(
  const dealii::Mapping<dim>        &mapping,
  const dealii::DoFHandler<dim>     &dof_handler,
  dealii::AffineConstraints<number> &constraints,
  boundaryCondition::type            boundary_type,
  unsigned int                       boundary_id,
  unsigned int                       component) const
{
  constexpr bool is_vector_field = spacedim != 1;

  // Create a component mask. This will only select a certain component for vector
  // fields
  dealii::ComponentMask mask = {};
  if constexpr (is_vector_field)
    {
      std::vector<bool> temp_mask(dim, false);
      temp_mask.at(component) = true;

      mask = dealii::ComponentMask(temp_mask);
    }

  // Apply the boundary conditions
  if (boundary_type == boundaryCondition::type::NATURAL)
    {
      // Do nothing because they are naturally enforced.
      return;
    }
  if (boundary_type == boundaryCondition::type::DIRICHLET)
    {
      dealii::VectorTools::interpolate_boundary_values(
        mapping,
        dof_handler,
        boundary_id,
        dealii::Functions::ZeroFunction<dim, number>(is_vector_field ? dim : 1),
        constraints,
        mask);
      return;
    }
  if (boundary_type == boundaryCondition::type::PERIODIC)
    {
      // Skip boundary ids that are odd since those map to the even faces
      if (boundary_id % 2 != 0)
        {
          return;
        }
      // Create a vector of matched pairs that we fill and enforce upon the
      // constaints
      std::vector<dealii::GridTools::PeriodicFacePair<
        typename dealii::DoFHandler<dim>::cell_iterator>>
        periodicity_vector;

      // Determine the direction
      const auto direction = static_cast<unsigned int>(std::floor(boundary_id / dim));

      // Collect the matched pairs on the coarsest level of the mesh
      dealii::GridTools::collect_periodic_faces(dof_handler,
                                                boundary_id,
                                                boundary_id + 1,
                                                direction,
                                                periodicity_vector);

      // Set constraints
      dealii::DoFTools::make_periodicity_constraints<dim, dim, number>(periodicity_vector,
                                                                       constraints,
                                                                       mask);
      return;
    }
  if (boundary_type == boundaryCondition::type::NEUMANN)
    {
      Assert(false, FeatureNotImplemented("Neumann boundary conditions"));
      return;
    }
  if (boundary_type == boundaryCondition::type::NON_UNIFORM_DIRICHLET)
    {
      dealii::VectorTools::interpolate_boundary_values(
        mapping,
        dof_handler,
        boundary_id,
        dealii::Functions::ZeroFunction<dim, number>(is_vector_field ? dim : 1),
        constraints,
        mask);
      return;
    }
  if (boundary_type == boundaryCondition::type::NON_UNIFORM_NEUMANN)
    {
      Assert(false, FeatureNotImplemented("Nonuniform neumann boundary conditions"));
      return;
    }
}

INSTANTIATE_UNI_TEMPLATE(constraintHandler)

PRISMS_PF_END_NAMESPACE
