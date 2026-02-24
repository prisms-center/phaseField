// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/numerics/vector_tools_boundary.h>

#include <prismspf/core/constraint_manager.h>
#include <prismspf/core/dof_manager.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/nonuniform_dirichlet.h>
#include <prismspf/core/system_wide.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/constraint_parameters.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
ConstraintManager<dim, degree, number>::ConstraintManager(
  const std::vector<FieldAttributes>         &field_attributes,
  const DofManager<dim>                      &_dof_manager,
  const PDEOperatorBase<dim, degree, number> *_pde_operator)
  : dof_manager(&_dof_manager)
  , pde_operator(_pde_operator)
  , constraints(field_attributes.size())
{
  /* TODO (fractalsbyx) figure out mg depth. Careful. Aux fields inherit this from
   * primary fields, as well as any dependencies*/
  unsigned int num_mg_levels = 1;
  for (unsigned int field_index = 0; field_index < field_attributes.size(); ++field_index)
    {
      constraints[field_index].resize(num_mg_levels);
      for (unsigned int relative_level = 0; relative_level < num_mg_levels;
           relative_level++)
        {
          const dealii::DoFHandler<dim> &dof_handler =
            dof_manager->get_dof_handler(field_index, relative_level);
          constraints[field_index][relative_level] = dealii::AffineConstraints<number>(
            dof_handler.locally_owned_dofs(),
            dealii::DoFTools::extract_locally_relevant_dofs(dof_handler));
        }
      // TODO (fractalsbyx): construct change_constraints
    }
}

template <unsigned int dim, unsigned int degree, typename number>
std::vector<const dealii::AffineConstraints<number> *>
ConstraintManager<dim, degree, number>::get_constraints(
  const std::set<Types::Index> &field_indices,
  unsigned int                  relative_level) const
{
  std::vector<const dealii::AffineConstraints<number> *> selected_constraints;
  selected_constraints.reserve(field_indices.size());
  for (const auto index : field_indices)
    {
      selected_constraints.push_back(&constraints[index][relative_level]);
    }
  return selected_constraints;
}

template <unsigned int dim, unsigned int degree, typename number>
const dealii::AffineConstraints<number> &
ConstraintManager<dim, degree, number>::get_constraint(Types::Index index,
                                                       unsigned int relative_level) const
{
  return constraints[index][relative_level];
}

template <unsigned int dim, unsigned int degree, typename number>
std::vector<const dealii::AffineConstraints<number> *>
ConstraintManager<dim, degree, number>::get_change_constraints(
  const std::set<Types::Index> &field_indices,
  unsigned int                  relative_level) const
{
  std::vector<const dealii::AffineConstraints<number> *> selected_constraints;
  selected_constraints.reserve(field_indices.size());
  for (const auto index : field_indices)
    {
      selected_constraints.push_back(&change_constraints[index][relative_level]);
    }
  return selected_constraints;
}

template <unsigned int dim, unsigned int degree, typename number>
const dealii::AffineConstraints<number> &
ConstraintManager<dim, degree, number>::get_change_constraint(
  Types::Index index,
  unsigned int relative_level) const
{
  return change_constraints[index][relative_level];
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::reinit(
  const std::vector<FieldAttributes> &field_attributes)
{
  const BoundaryParameters<dim> &boundary_parameters =
    user_inputs->get_boundary_parameters();
  // The map from user inputs has string keys for now.
  std::map<std::string, Types::Index> field_indices = field_index_map(field_attributes);
  for (const auto &[name, field_constraints] :
       boundary_parameters.boundary_condition_list)
    {
      const unsigned int field_index = field_indices.at(name);
      for (unsigned int relative_level = 0;
           relative_level < constraints[field_index].size();
           ++relative_level)
        {
          dealii::AffineConstraints<number> &constraint =
            constraints[field_index][relative_level];
          dealii::AffineConstraints<number> &change_constraint =
            change_constraints[field_index][relative_level];
          const dealii::DoFHandler<dim> &dof_handler =
            dof_manager->get_dof_handler(field_index, relative_level);

          make_constraints_for_single_field(constraint,
                                            dof_handler,
                                            field_constraints,
                                            field_attributes[field_index].field_type,
                                            field_index,
                                            false);
          make_constraints_for_single_field(change_constraint,
                                            dof_handler,
                                            field_constraints,
                                            field_attributes[field_index].field_type,
                                            field_index,
                                            true);
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::make_constraints_for_single_field(
  dealii::AffineConstraints<number> &constraint,
  const dealii::DoFHandler<dim>     &dof_handler,
  const FieldConstraints<dim>       &field_constraints,
  TensorRank                         tensor_rank,
  Types::Index                       field_index,
  bool                               for_change_term)
{
  // 0. Reinitialize constraint with the correct dof numbering
  constraint.reinit(dof_handler.locally_owned_dofs(),
                    dealii::DoFTools::extract_locally_relevant_dofs(dof_handler));

  // 1. Make hanging node constraints
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraint);

  // 2. Make boundary constraints
  make_bc_constraints(constraint,
                      dof_handler,
                      field_constraints,
                      tensor_rank,
                      field_index,
                      for_change_term);

  // 3. TODO: Make pinned point constraints
  // if (boundary_parameters.has_pinned_point(field_index))
  //  {
  //    const auto &[value, target_point] =
  //      boundary_parameters.get_pinned_point(field_index);
  //    set_pinned_point(constraint, target_point, value, dof_handler, for_change_term);
  //  }
  constraint.close();
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::make_bc_constraints(
  dealii::AffineConstraints<number> &constraint,
  const dealii::DoFHandler<dim>     &dof_handler,
  const FieldConstraints<dim>       &boundary_condition,
  const TensorRank                   tensor_rank,
  Types::Index                       field_index,
  bool                               for_change_term)
{
  for (unsigned int comp = 0; comp < dim; comp++)
    {
      const ComponentConditions &comp_bcs =
        boundary_condition.component_constraints.at(comp);
      for (const auto &[boundary_id, boundary_type] : comp_bcs.conditions)
        {
          make_one_boundary_constraint(constraint,
                                       boundary_id,
                                       comp,
                                       boundary_type,
                                       dof_handler,
                                       tensor_rank,
                                       field_index,
                                       for_change_term);
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::make_one_boundary_constraint(
  dealii::AffineConstraints<number> &_constraints,
  unsigned int                       boundary_id,
  unsigned int                       component,
  Condition                          boundary_type,
  const dealii::DoFHandler<dim>     &dof_handler,
  TensorRank                         tensor_rank,
  Types::Index                       field_index,
  bool                               for_change_term) const
{
  const bool                  is_vector_field = tensor_rank == TensorRank::Vector;
  const dealii::ComponentMask mask =
    is_vector_field ? vector_component_mask.at(component) : scalar_empty_mask;

  // Apply the boundary conditions
  switch (boundary_type)
    {
      case Condition::Natural:
        {
          // do nothing
          break;
        }
      case Condition::Dirichlet:
      case Condition::TimeDependentDirichlet:
        {
          make_dirichlet_constraints(_constraints,
                                     dof_handler,
                                     boundary_id,
                                     field_index,
                                     is_vector_field,
                                     mask,
                                     for_change_term);
          break;
        }
      case Condition::UniformDirichlet:
        {
          Assert(false, FeatureNotImplemented("Uniform dirichlet boundary conditions"));
          break;
        }
      case Condition::Neumann:
        {
          Assert(false, FeatureNotImplemented("Neumann boundary conditions"));
          break;
        }

      case Condition::TimeDependentNeumann:
        {
          Assert(false,
                 FeatureNotImplemented("Time dependent neumann boundary conditions"));
          break;
        }
      case Condition::UniformNeumann:
        {
          Assert(false, FeatureNotImplemented("Uniform neumann boundary conditions"));
          break;
        }
      default:
        {
          AssertThrow(false, UnreachableCode());
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::update_time_dependent_constraints(
  const std::vector<FieldAttributes> &field_attributes)
{
  const BoundaryParameters<dim> &boundary_parameters =
    user_inputs->get_boundary_parameters();
  // The map from user inputs has string keys for now.
  std::map<std::string, Types::Index> field_indices = field_index_map(field_attributes);
  for (const auto &[name, field_constraints] :
       boundary_parameters.boundary_condition_list)
    {
      if (field_constraints.has_time_dependent_bcs())
        {
          const unsigned int field_index = field_indices.at(name);
          for (unsigned int relative_level = 0;
               relative_level < constraints[field_index].size();
               ++relative_level)
            {
              dealii::AffineConstraints<number> &constraint =
                constraints[field_index][relative_level];
              dealii::AffineConstraints<number> &change_constraint =
                change_constraints[field_index][relative_level];
              const dealii::DoFHandler<dim> &dof_handler =
                dof_manager->get_dof_handler(field_index, relative_level);

              make_constraints_for_single_field(constraint,
                                                dof_handler,
                                                field_constraints,
                                                field_attributes[field_index].field_type,
                                                field_index,
                                                false);
              make_constraints_for_single_field(change_constraint,
                                                dof_handler,
                                                field_constraints,
                                                field_attributes[field_index].field_type,
                                                field_index,
                                                true);
            }
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
const std::array<dealii::ComponentMask, dim>
  ConstraintManager<dim, degree, number>::vector_component_mask = []()
{
  std::array<dealii::ComponentMask, dim> masks {};
  for (unsigned int i = 0; i < dim; ++i)
    {
      dealii::ComponentMask temp_mask(dim, false);
      temp_mask.set(i, true);
      masks.at(i) = temp_mask;
    }
  return masks;
}();

template <unsigned int dim, unsigned int degree, typename number>
const dealii::ComponentMask ConstraintManager<dim, degree, number>::scalar_empty_mask {};

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::make_natural_constraints() const
{
  // Do nothing because they are naturally enforced.
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConstraintManager<dim, degree, number>::make_dirichlet_constraints(
  dealii::AffineConstraints<number> &_constraints,
  const dealii::DoFHandler<dim>     &dof_handler,
  const unsigned int                &boundary_id,
  const unsigned int                &field_index,
  const bool                        &is_vector_field,
  const dealii::ComponentMask       &mask,
  bool                               is_change_term) const
{
  if (!is_change_term)
    {
      dealii::VectorTools::interpolate_boundary_values(
        SystemWide<dim, degree>::mapping,
        dof_handler,
        boundary_id,
        NonuniformDirichlet<dim, degree, number>(field_index,
                                                 boundary_id,
                                                 *pde_operator,
                                                 is_vector_field ? dim : 1),
        _constraints,
        mask);
    }
  else
    {
      dealii::VectorTools::interpolate_boundary_values(
        SystemWide<dim, degree>::mapping,
        dof_handler,
        boundary_id,
        dealii::Functions::ZeroFunction<dim, number>(is_vector_field ? dim : 1),
        _constraints,
        mask);
    }
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
  dealii::AffineConstraints<number> &constraint,
  const dealii::Point<dim>          &target_point,
  const std::array<number, dim>     &value,
  const dealii::DoFHandler<dim>     &dof_handler,
  const TensorRank                   tensor_rank,
  bool                               is_change_term) const
{
  const double tolerance = 1.0e-2;

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
          const unsigned int dimension = tensor_rank == TensorRank::Vector ? dim : 1;
          for (unsigned int component = 0; component < dimension; ++component)
            {
              constraint.add_line(dof_index + component);
              constraint.set_inhomogeneity(dof_index + component,
                                           is_change_term ? value.at(component) : 0.0);
            }
        }
    }
}

#include "core/constraint_manager.inst"

PRISMS_PF_END_NAMESPACE
