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

#include <prismspf/core/constraint_manager.h>
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
#include <string>
#include <type_traits>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
ConstraintManager<dim, degree, number>::ConstraintManager(
  const std::vector<FieldAttributes>                            &field_attributes,
  const std::set<SolveGroup>                                    &solve_groups,
  const std::shared_ptr<const PDEOperator<dim, degree, number>> &_pde_operator)
  : pde_operator(_pde_operator)
  , constraints(field_attributes.size())
{
  for (const auto &solve_group : solve_groups)
    {
      /* TODO (fractalsbyx) figure out mg depth. Careful. Aux fields inherit this from
       * primary fields, as well as any dependencies*/
      unsigned int num_mg_levels = 1;
      for (const auto &field_index : solve_group.field_indices)
        {
          constraints[field_index].resize(num_mg_levels, nullptr);
          for (unsigned int relative_level = 0; relative_level < num_mg_levels;
               relative_level++)
            {
              constraints[field_index][relative_level] =
                std::make_shared<dealii::AffineConstraints<number>>();
            }
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
std::vector<const dealii::AffineConstraints<number> *>
ConstraintManager<dim, degree, number>::get_constraints(
  const std::set<unsigned int> &field_indices,
  unsigned int                  relative_level) const
{
  std::vector<const dealii::AffineConstraints<number> *> selected_constraints;
  selected_constraints.reserve(field_indices.size());
  for (const auto index : field_indices)
    {
      selected_constraints.push_back(constraints[index][relative_level].get());
    }
  return selected_constraints;
}

template <unsigned int dim, unsigned int degree, typename number>
const dealii::AffineConstraints<number> &
ConstraintManager<dim, degree, number>::get_constraint(Types::Index index,
                                                       unsigned int relative_level) const
{
  return *constraints[index][relative_level];
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::make_constraint(
  const dealii::Mapping<dim>    &mapping,
  const dealii::DoFHandler<dim> &dof_handler,
  unsigned int                   index)
{
  // The constraint that we are going to modify
  auto &local_constraint = constraints[index];

  apply_generic_constraints(dof_handler, local_constraint);

  // First check the normal boundary conditions
  const auto &boundary_condition =
    user_inputs->get_boundary_parameters().get_boundary_condition_list().at(index);
  for (const auto &[component, condition] : boundary_condition)
    {
      for (const auto &[boundary_id, boundary_type] :
           condition.get_boundary_condition_map())
        {
          if (user_inputs->get_variable_attributes().at(index).field_info.tensor_rank !=
              FieldInfo::TensorRank::Vector)
            {
              apply_constraints<number, 1>(mapping,
                                           dof_handler,
                                           local_constraint,
                                           condition,
                                           boundary_type,
                                           boundary_id,
                                           component,
                                           index);
              continue;
            }
          apply_constraints<number, dim>(mapping,
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
  if (user_inputs->get_boundary_parameters().has_pinned_point(index))
    {
      set_pinned_point<number>(dof_handler, local_constraint, index, false);
    }

  // Close constraints
  local_constraint.close();
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::make_constraints(
  const dealii::Mapping<dim>                         &mapping,
  const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers)
{
  for (const auto &[index, variable] : user_inputs->get_variable_attributes())
    {
      make_constraint(mapping, *dof_handlers.at(index), index);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::update_time_dependent_constraints(
  const dealii::Mapping<dim>                         &mapping,
  const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers)
{
  for (const auto &[index, variable] : user_inputs->get_variable_attributes())
    {
      // Check that we have time-dependent constraints before recreating the whole
      // constraint set.
      if (user_inputs->get_boundary_parameters().is_time_dependent(index))
        {
          // TODO (landinjm): Is there a way to update the constraint set without
          // recreating
          // it?
          make_constraint(mapping, *dof_handlers.at(index), index);
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::ComponentMask
ConstraintManager<dim, degree, number>::create_component_mask(unsigned int component,
                                                              bool is_vector_field) const
{
  if (!is_vector_field)
    {
      return {};
    }

  dealii::ComponentMask temp_mask(dim, false);
  temp_mask.set(component, true);
  return temp_mask;
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::make_natural_constraints() const
{
  // Do nothing because they are naturally enforced.
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::make_uniform_dirichlet_constraints(
  const dealii::Mapping<dim>        &mapping,
  const dealii::DoFHandler<dim>     &dof_handler,
  const unsigned int                &boundary_id,
  const bool                        &is_vector_field,
  const number                      &value,
  dealii::AffineConstraints<number> &_constraints,
  const dealii::ComponentMask       &mask) const
{
  dealii::VectorTools::interpolate_boundary_values(
    mapping,
    dof_handler,
    boundary_id,
    dealii::Functions::ConstantFunction<dim, number>(value, is_vector_field ? dim : 1),
    _constraints,
    mask);
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::make_nonuniform_dirichlet_constraints(
  const dealii::Mapping<dim>        &mapping,
  const dealii::DoFHandler<dim>     &dof_handler,
  const unsigned int                &boundary_id,
  const unsigned int                &index,
  const bool                        &is_vector_field,
  dealii::AffineConstraints<number> &_constraints,
  const dealii::ComponentMask       &mask,
  bool                               is_change_term) const
{
  if (!is_change_term)
    {
      dealii::VectorTools::interpolate_boundary_values(
        mapping,
        dof_handler,
        boundary_id,
        NonuniformDirichlet<dim, degree, number>(index,
                                                 boundary_id,
                                                 pde_operator,
                                                 is_vector_field ? dim : 1),
        _constraints,
        mask);
      return;
    }

  dealii::VectorTools::interpolate_boundary_values(
    mapping,
    dof_handler,
    boundary_id,
    dealii::Functions::ZeroFunction<dim, number>(is_vector_field ? dim : 1),
    _constraints,
    mask);
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::make_periodic_constraints(
  dealii::AffineConstraints<number> &_constraints,
  const dealii::DoFHandler<dim>     &dof_handler,
  const unsigned int                &boundary_id,
  const dealii::ComponentMask       &mask) const
{
  // Create a vector of matched pairs that we fill and enforce upon the
  // constaints
  std::vector<
    dealii::GridTools::PeriodicFacePair<typename dealii::DoFHandler<dim>::cell_iterator>>
    periodicity_vector;

  // Determine the direction
  const unsigned int direction = boundary_id / 2;

  // Collect the matched pairs on the coarsest level of the mesh
  dealii::GridTools::collect_periodic_faces(dof_handler,
                                            boundary_id,
                                            boundary_id + 1,
                                            direction,
                                            periodicity_vector);

  // Set constraints
  dealii::DoFTools::make_periodicity_constraints<dim, dim, number>(periodicity_vector,
                                                                   _constraints,
                                                                   mask);
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::set_pinned_point(
  const dealii::DoFHandler<dim>     &dof_handler,
  dealii::AffineConstraints<number> &_constraints,
  const dealii::Point<dim>          &target_point,
  unsigned int                       index,
  bool                               is_change_term) const
{
  const double tolerance = 1.0e-2;
  const auto &[value, target_point] =
    user_inputs->get_boundary_parameters().get_pinned_point(index);

  // Helper function to set inhomogeneity for a single DOF
  auto set_inhomogeneity = [&](unsigned int dof_index, double _value)
  {
    _constraints.add_line(dof_index);
    _constraints.set_inhomogeneity(dof_index, is_change_term ? _value : 0.0);
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

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::apply_generic_constraints(
  const dealii::DoFHandler<dim>     &dof_handler,
  dealii::AffineConstraints<number> &_constraints) const
{
  // Clear constraints
  _constraints.clear();

  // Reinitialize constraints
  _constraints.reinit(dof_handler.locally_owned_dofs(),
                      dealii::DoFTools::extract_locally_relevant_dofs(dof_handler));

  // Make hanging node constraints
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, _constraints);
}

template <unsigned int dim, unsigned int degree, typename number>
template <int spacedim>
void
ConstraintManager<dim, degree, number>::apply_constraints(
  const unsigned int                 index,
  const unsigned int                 boundary_id,
  const unsigned int                 component,
  dealii::AffineConstraints<number> &_constraints,
  BoundaryCondition::Type            boundary_type,
  const number                       dirichlet_value,
  const dealii::Mapping<dim>        &mapping,
  const dealii::DoFHandler<dim>     &dof_handler,
  bool                               for_change_term) const
{
  constexpr bool              is_vector_field = spacedim != 1;
  const dealii::ComponentMask mask = create_component_mask(component, is_vector_field);

  // Apply the boundary conditions
  switch (boundary_type)
    {
      case BoundaryCondition::Type::Natural:
        {
          make_natural_constraints();
          break;
        }
      case BoundaryCondition::Type::Dirichlet:
        {
          make_uniform_dirichlet_constraints(mapping,
                                             dof_handler,
                                             boundary_id,
                                             is_vector_field,
                                             for_change_term ? 0.0 : dirichlet_value,
                                             _constraints,
                                             mask);
          break;
        }
      case BoundaryCondition::Type::NonuniformDirichlet:
      case BoundaryCondition::Type::TimeDependentNonuniformDirichlet:
        {
          make_nonuniform_dirichlet_constraints(mapping,
                                                dof_handler,
                                                boundary_id,
                                                index,
                                                is_vector_field,
                                                _constraints,
                                                mask,
                                                for_change_term);
          break;
        }
      case BoundaryCondition::Type::Periodic:
        {
          // Skip boundary ids that are odd since those map to the even faces
          if (boundary_id % 2 != 0)
            {
              break;
            }
          make_periodic_constraints(dof_handler, boundary_id, _constraints, mask);
          break;
        }
      case BoundaryCondition::Type::Neumann:
        {
          Assert(false, FeatureNotImplemented("Neumann boundary conditions"));
          break;
        }
      case BoundaryCondition::Type::NonuniformNeumann:
        {
          Assert(false, FeatureNotImplemented("Nonuniform neumann boundary conditions"));
          break;
        }
      case BoundaryCondition::Type::TimeDependentNonuniformNeumann:
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

#include "core/constraint_handler.inst"

PRISMS_PF_END_NAMESPACE
