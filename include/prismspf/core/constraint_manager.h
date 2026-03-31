// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/affine_constraints.h>

#include <prismspf/core/dof_manager.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_block.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/constraint_parameters.h>
#include <prismspf/user_inputs/spatial_discretization.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

// Forward declaration to avoid circular dependency.
template <unsigned int dim, unsigned int degree, typename number>
class PDEOperatorBase;

// TODO (fractalsbyx): The following snippet is from dealii. Using this (by
// pre-constructing the needed maps and functions) may be cleaner than how dirichlet are
// currently done.
/* template <int dim, int spacedim, typename number>
   void
  interpolate_boundary_values(
    const Mapping<dim, spacedim>    &mapping,
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                              &function_map,
    AffineConstraints<number> &constraints,
    const ComponentMask       &component_mask = {}); */

// TODO (fractalsbyx): I don't think we actually need 'change' constraints for newton
// solves at all. Verify this and remove them.

/**
 * @brief The class handles the generation and application of boundary conditions based on
 * the user-inputs.
 */
template <unsigned int dim, unsigned int degree, typename number>
class ConstraintManager
{
public:
  /**
   * @brief Constructor.
   */
  ConstraintManager(const std::vector<FieldAttributes>         &field_attributes,
                    const BoundaryParameters<dim>              &_boundary_parameters,
                    const SpatialDiscretization<dim>           &_spatial_discretization,
                    const DoFManager<dim, degree>              &_dof_manager,
                    const PDEOperatorBase<dim, degree, number> &_pde_operator);

  /**
   * @brief Getter function for the constraints.
   */
  [[nodiscard]] std::vector<const dealii::AffineConstraints<number> *>
  get_constraints(const std::set<unsigned int> &field_indices,
                  unsigned int                  relative_level = 0) const;

  /**
   * @brief Getter function for the constraint of an index (constant reference).
   */
  [[nodiscard]] const dealii::AffineConstraints<number> &
  get_constraint(Types::Index index, unsigned int relative_level = 0) const;

  /**
   * @brief Getter function for the constraints.
   */
  [[nodiscard]] const std::vector<std::array<dealii::AffineConstraints<number>, 2>> &
  get_generic_constraints() const;

  /**
   * @brief Getter function for the constraint of an index (constant reference).
   */
  [[nodiscard]] const dealii::AffineConstraints<number> &
  get_generic_constraint(unsigned int rank, unsigned int relative_level = 0) const;

  /**
   * @brief Make constraints based on the inputs of the constructor.
   */
  void
  reinit(const std::vector<FieldAttributes> &field_attributes);

  /**
   * @brief Update time-dependent constraints.
   * For now this only updates the non-uniform dirichlet constraints.
   */
  void
  update_time_dependent_constraints(const std::vector<FieldAttributes> &field_attributes);

private:
  /**
   * @brief Create a component mask.
   */
  static const std::array<dealii::ComponentMask, dim> vector_component_mask;
  static const dealii::ComponentMask                  scalar_empty_mask;

  /**
   * @brief Construct constraints for a single field based on the boundary conditions.
   */
  void
  make_constraints_for_single_field(dealii::AffineConstraints<number> &constraint,
                                    const dealii::DoFHandler<dim>     &dof_handler,
                                    const FieldConstraints<dim>       &field_constraints,
                                    TensorRank                         tensor_rank,
                                    Types::Index                       field_index);

  /**
   * @brief Add boundary conditions to a single constraint.
   */
  void
  make_bc_constraints(dealii::AffineConstraints<number> &constraint,
                      const dealii::DoFHandler<dim>     &dof_handler,
                      const FieldConstraints<dim>       &boundary_condition,
                      TensorRank                         tensor_rank,
                      Types::Index                       field_index);

  /**
   * @brief Apply constraints for common boundary conditions.
   */
  void
  make_one_boundary_constraint(dealii::AffineConstraints<number> &_constraints,
                               unsigned int                       boundary_id,
                               unsigned int                       component,
                               Condition                          boundary_type,
                               const dealii::DoFHandler<dim>     &dof_handler,
                               TensorRank                         tensor_rank,
                               Types::Index                       field_index) const;

  /**
   * @brief Apply natural constraints.
   */
  void
  make_natural_constraints() const;

  /**
   * @brief Make dirichlet constraints.
   */
  void
  make_dirichlet_constraints(dealii::AffineConstraints<number> &_constraints,
                             const dealii::DoFHandler<dim>     &dof_handler,
                             const unsigned int                &boundary_id,
                             const unsigned int                &field_index,
                             const bool                        &is_vector_field,
                             const dealii::ComponentMask       &mask) const;

  /**
   * @deprecated We apply periodic conditions to the mesh itself
   * @brief make periodic constraints.
   */
  void
  make_periodic_constraints(dealii::AffineConstraints<number> &_constraints,
                            const dealii::DoFHandler<dim>     &dof_handler,
                            const unsigned int                &boundary_id,
                            const dealii::ComponentMask       &mask) const;

  /**
   * @brief Set the dirichlet constraint for the pinned point.
   */
  void
  set_pinned_point(dealii::AffineConstraints<number> &constraint,
                   const dealii::Point<dim>          &target_point,
                   const std::array<number, dim>     &value,
                   const dealii::DoFHandler<dim>     &dof_handler,
                   TensorRank                         tensor_rank) const;

  /**
   * @brief User-inputs constraint parameters.
   */
  const BoundaryParameters<dim> *boundary_parameters;

  /**
   * @brief User-inputs discretization.
   */
  const SpatialDiscretization<dim> *spatial_discretization;

  /**
   * @brief Dof manager pointer.
   */
  const DoFManager<dim, degree> *dof_manager;

  /**
   * @brief PDE operator.
   */
  const PDEOperatorBase<dim, degree, number> *pde_operator;

  /**
   * @brief Constraints. Outer vector is indexed by field index. Inner vector is indexed
   * by relative mg level.
   */
  std::vector<std::vector<dealii::AffineConstraints<number>>> constraints;

  /**
   * @brief Constraints not specific to any field. We need this for invm
   */
  std::vector<std::array<dealii::AffineConstraints<number>, 2>> generic_constraints;
};

PRISMS_PF_END_NAMESPACE
