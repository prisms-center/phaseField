// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/numerics/vector_tools_boundary.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/nonuniform_dirichlet.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/boundary_parameters.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
constraintHandler<dim, degree>::constraintHandler(
  const userInputParameters<dim>                                &_user_inputs,
  const MGInfo<dim>                                             &_mg_info,
  const std::shared_ptr<const PDEOperator<dim, degree, double>> &_pde_operator,
  const std::shared_ptr<const PDEOperator<dim, degree, float>>  &_pde_operator_float)
  : user_inputs(&_user_inputs)
  , mg_info(&_mg_info)
  , pde_operator(_pde_operator)
  , pde_operator_float(_pde_operator_float)
  , constraints(_user_inputs.get_variable_attributes().size(),
                dealii::AffineConstraints<double>())
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

template <unsigned int dim, unsigned int degree>
std::vector<const dealii::AffineConstraints<double> *>
constraintHandler<dim, degree>::get_constraints()
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

template <unsigned int dim, unsigned int degree>
const dealii::AffineConstraints<double> &
constraintHandler<dim, degree>::get_constraint(unsigned int index) const
{
  Assert(constraints.size() > index,
         dealii::ExcMessage("The constraint set does not contain index = " +
                            std::to_string(index)));
  return constraints.at(index);
}

template <unsigned int dim, unsigned int degree>
std::vector<const dealii::AffineConstraints<float> *>
constraintHandler<dim, degree>::get_mg_constraints(unsigned int level)
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

template <unsigned int dim, unsigned int degree>
const dealii::AffineConstraints<float> &
constraintHandler<dim, degree>::get_mg_constraint(unsigned int level,
                                                  unsigned int index) const
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

template <unsigned int dim, unsigned int degree>
void
constraintHandler<dim, degree>::make_constraints(
  const dealii::Mapping<dim>                         &mapping,
  const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers)
{
  for (const auto &[index, variable] : user_inputs->get_variable_attributes())
    {
      make_constraint(mapping, *dof_handlers.at(index), index);
    }
}

template <unsigned int dim, unsigned int degree>
void
constraintHandler<dim, degree>::make_mg_constraints(
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

template <unsigned int dim, unsigned int degree>
void
constraintHandler<dim, degree>::update_time_dependent_constraints(
  const dealii::Mapping<dim>                         &mapping,
  const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers)
{
  for (const auto &[index, variable] : user_inputs->get_variable_attributes())
    {
      // Check that we have time-dependent constraints before recreating the whole
      // constraint set.
      if (user_inputs->get_boundary_parameters().time_dependent_BC_list.contains(index))
        {
          // TODO (landinjm): Is there a way to update the constraint set without
          // recreating
          // it?
          make_constraint(mapping, *dof_handlers.at(index), index);
        }
    }
}

template <unsigned int dim, unsigned int degree>
void
constraintHandler<dim, degree>::update_time_dependent_mg_constraints(
  const dealii::Mapping<dim>                         &mapping,
  const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers,
  unsigned int                                        level)
{
  Assert(has_multigrid, dealii::ExcNotInitialized());

  for (unsigned int index = 0; index < dof_handlers.size(); index++)
    {
      if (user_inputs->get_boundary_parameters().time_dependent_BC_list.contains(index))
        {
          // TODO (landinjm): Fix this so we can actually apply constraints to the LHS
          // fields baseded on whether it is normal or change. For now, I know this will
          // have no effect on the precipitate app, so I'll leave it.

          // TODO (landinjm): Is there a way to update the constraint set without
          // recreating it?
          make_mg_constraint(mapping,
                             *dof_handlers[index],
                             index,
                             level,
                             dependencyType::CHANGE);
        }
    }
}

template <unsigned int dim, unsigned int degree>
dealii::ComponentMask
constraintHandler<dim, degree>::create_component_mask(unsigned int component,
                                                      bool         is_vector_field) const
{
  if (!is_vector_field)
    {
      return {};
    }

  dealii::ComponentMask temp_mask(dim, false);
  temp_mask.set(component, true);
  return temp_mask;
}

template <unsigned int dim, unsigned int degree>
void
constraintHandler<dim, degree>::apply_natural_constraints() const
{
  // Do nothing because they are naturally enforced.
}

template <unsigned int dim, unsigned int degree>
template <typename number>
void
constraintHandler<dim, degree>::apply_dirichlet_constraints(
  const dealii::Mapping<dim>        &mapping,
  const dealii::DoFHandler<dim>     &dof_handler,
  const unsigned int                &boundary_id,
  const bool                        &is_vector_field,
  const number                      &value,
  dealii::AffineConstraints<number> &constraints,
  const dealii::ComponentMask       &mask) const
{
  dealii::VectorTools::interpolate_boundary_values(
    mapping,
    dof_handler,
    boundary_id,
    dealii::Functions::ConstantFunction<dim, number>(value, is_vector_field ? dim : 1),
    constraints,
    mask);
}

template <unsigned int dim, unsigned int degree>
template <typename number>
void
constraintHandler<dim, degree>::apply_periodic_constraints(
  const dealii::DoFHandler<dim>     &dof_handler,
  const unsigned int                &boundary_id,
  dealii::AffineConstraints<number> &constraints,
  const dealii::ComponentMask       &mask) const
{
  // Create a vector of matched pairs that we fill and enforce upon the
  // constaints
  std::vector<
    dealii::GridTools::PeriodicFacePair<typename dealii::DoFHandler<dim>::cell_iterator>>
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
}

template <unsigned int dim, unsigned int degree>
template <typename number>
void
constraintHandler<dim, degree>::apply_nonuniform_dirichlet_constraints(
  const dealii::Mapping<dim>        &mapping,
  const dealii::DoFHandler<dim>     &dof_handler,
  const unsigned int                &boundary_id,
  const unsigned int                &index,
  const bool                        &is_vector_field,
  dealii::AffineConstraints<number> &constraints,
  const dealii::ComponentMask       &mask,
  bool                               is_change_term) const
{
  if (!is_change_term)
    {
      const auto &desired_pde_operator = [this]() -> const auto &
      {
        if constexpr (std::is_same_v<number, double>)
          {
            return pde_operator;
          }
        else
          {
            return pde_operator_float;
          }
      }();

      dealii::VectorTools::interpolate_boundary_values(
        mapping,
        dof_handler,
        boundary_id,
        nonuniformDirichlet<dim, degree, number>(index,
                                                 boundary_id,
                                                 desired_pde_operator,
                                                 is_vector_field ? dim : 1),
        constraints,
        mask);
      return;
    }

  dealii::VectorTools::interpolate_boundary_values(
    mapping,
    dof_handler,
    boundary_id,
    dealii::Functions::ZeroFunction<dim, number>(is_vector_field ? dim : 1),
    constraints,
    mask);
}

template <unsigned int dim, unsigned int degree>
void
constraintHandler<dim, degree>::make_constraint(
  const dealii::Mapping<dim>    &mapping,
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
    user_inputs->get_boundary_parameters().boundary_condition_list.at(index);
  for (const auto &[component, condition] : boundary_condition)
    {
      for (const auto &[boundary_id, boundary_type] :
           condition.get_boundary_condition_map())
        {
          if (user_inputs->get_variable_attributes().at(index).field_type !=
              fieldType::VECTOR)
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
  if (user_inputs->get_boundary_parameters().pinned_point_list.contains(index))
    {
      set_pinned_point<double>(dof_handler, local_constraint, index, false);
    }

  // Close constraints
  local_constraint.close();
}

template <unsigned int dim, unsigned int degree>
void
constraintHandler<dim, degree>::make_mg_constraint(
  const dealii::Mapping<dim>    &mapping,
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
        this->user_inputs->get_boundary_parameters().boundary_condition_list.at(
          global_index);
      for (const auto &[component, condition] : boundary_condition)
        {
          for (const auto &[boundary_id, boundary_type] :
               condition.get_boundary_condition_map())
            {
              if (user_inputs->get_variable_attributes().at(global_index).field_type !=
                  fieldType::VECTOR)
                {
                  apply_constraints<float, 1>(mapping,
                                              dof_handler,
                                              local_constraint,
                                              condition,
                                              boundary_type,
                                              boundary_id,
                                              component,
                                              global_index,
                                              true);
                  continue;
                }
              apply_constraints<float, dim>(mapping,
                                            dof_handler,
                                            local_constraint,
                                            condition,
                                            boundary_type,
                                            boundary_id,
                                            component,
                                            global_index,
                                            true);
            }
        }

      // Second check for pinned points, if they exist
      if (user_inputs->get_boundary_parameters().pinned_point_list.contains(global_index))
        {
          set_pinned_point<float>(dof_handler, local_constraint, global_index, true);
        }
    }
  else if (dependency_type == dependencyType::NORMAL)
    {
      const auto &boundary_condition =
        this->user_inputs->get_boundary_parameters().boundary_condition_list.at(
          global_index);
      for (const auto &[component, condition] : boundary_condition)
        {
          for (const auto &[boundary_id, boundary_type] :
               condition.get_boundary_condition_map())
            {
              if (user_inputs->get_variable_attributes().at(global_index).field_type !=
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

      // Second check for pinned points, if they exist
      if (user_inputs->get_boundary_parameters().pinned_point_list.contains(global_index))
        {
          set_pinned_point<float>(dof_handler, local_constraint, global_index, false);
        }
    }
  else
    {
      AssertThrow(false, dealii::ExcNotImplemented());
    }

  local_constraint.close();
}

template <unsigned int dim, unsigned int degree>
template <typename number>
void
constraintHandler<dim, degree>::set_pinned_point(
  const dealii::DoFHandler<dim>     &dof_handler,
  dealii::AffineConstraints<number> &constraints,
  unsigned int                       index,
  bool                               is_change_term) const
{
  const number tolerance = 1.0e-2;
  const auto &[value, target_point] =
    user_inputs->get_boundary_parameters().pinned_point_list.at(index);

  // Helper function to set inhomogeneity for a single DOF
  auto set_inhomogeneity = [&](unsigned int dof_index, number value)
  {
    constraints.add_line(dof_index);
    constraints.set_inhomogeneity(dof_index, is_change_term ? value : 0.0);
  };

  // Helper function to handle vector values
  auto set_vector_inhomogeneity =
    [&](unsigned int base_dof_index, const std::vector<double> &values)
  {
    for (unsigned int dimension = 0; dimension < dim; ++dimension)
      {
        set_inhomogeneity(base_dof_index + dimension, values[dimension]);
      }
  };

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        {
          continue;
        }

      for (unsigned int vertex = 0; vertex < dealii::GeometryInfo<dim>::vertices_per_cell;
           ++vertex)
        {
          const auto vertex_point = cell->vertex(vertex);
          if (target_point.distance(vertex_point) >= tolerance * cell->diameter())
            {
              continue;
            }

          const unsigned int dof_index = cell->vertex_dof_index(vertex, 0);

          // Handle both scalar and vector values
          std::visit(
            [&](const auto &val)
            {
              if constexpr (std::is_same_v<std::decay_t<decltype(val)>, double>)
                {
                  set_inhomogeneity(dof_index, val);
                }
              else
                {
                  set_vector_inhomogeneity(dof_index, val);
                }
            },
            value);
        }
    }
}

template <unsigned int dim, unsigned int degree>
template <typename number>
void
constraintHandler<dim, degree>::apply_generic_constraints(
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

template <unsigned int dim, unsigned int degree>
template <typename number, int spacedim>
void
constraintHandler<dim, degree>::apply_constraints(
  const dealii::Mapping<dim>        &mapping,
  const dealii::DoFHandler<dim>     &dof_handler,
  dealii::AffineConstraints<number> &constraints,
  const boundaryCondition           &boundary_condition,
  boundaryCondition::type            boundary_type,
  unsigned int                       boundary_id,
  unsigned int                       component,
  unsigned int                       index,
  bool                               is_change_term) const
{
  constexpr bool              is_vector_field = spacedim != 1;
  const dealii::ComponentMask mask = create_component_mask(component, is_vector_field);

  // Apply the boundary conditions
  switch (boundary_type)
    {
      case boundaryCondition::type::NATURAL:
        {
          apply_natural_constraints();
          break;
        }
      case boundaryCondition::type::DIRICHLET:
        {
          apply_dirichlet_constraints<number>(
            mapping,
            dof_handler,
            boundary_id,
            is_vector_field,
            is_change_term ? 0.0 : boundary_condition.get_dirichlet_value(boundary_id),
            constraints,
            mask);
          break;
        }
      case boundaryCondition::type::PERIODIC:
        {
          // Skip boundary ids that are odd since those map to the even faces
          if (boundary_id % 2 != 0)
            {
              break;
            }
          apply_periodic_constraints<number>(dof_handler, boundary_id, constraints, mask);
          break;
        }
      case boundaryCondition::type::NEUMANN:
        {
          Assert(false, FeatureNotImplemented("Neumann boundary conditions"));
          break;
        }
      case boundaryCondition::type::NON_UNIFORM_DIRICHLET:
      case boundaryCondition::type::TIME_DEPENDENT_NON_UNIFORM_DIRICHLET:
        {
          apply_nonuniform_dirichlet_constraints<number>(mapping,
                                                         dof_handler,
                                                         boundary_id,
                                                         index,
                                                         is_vector_field,
                                                         constraints,
                                                         mask,
                                                         is_change_term);
          break;
        }
      case boundaryCondition::type::NON_UNIFORM_NEUMANN:
        {
          Assert(false, FeatureNotImplemented("Nonuniform neumann boundary conditions"));
          break;
        }
      case boundaryCondition::type::TIME_DEPENDENT_NON_UNIFORM_NEUMANN:
        {
          Assert(false,
                 FeatureNotImplemented(
                   "Time dependent nonuniform neumann boundary conditions"));
          break;
        }
      default:
        {
          AssertThrow(false, UnreachableCode());
        }
    }
}

INSTANTIATE_BI_TEMPLATE(constraintHandler)

PRISMS_PF_END_NAMESPACE
