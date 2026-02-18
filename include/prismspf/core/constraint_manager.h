// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/affine_constraints.h>

#include <prismspf/core/dof_manager.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/boundary_parameters.h>
#include <prismspf/user_inputs/user_input_parameters.h>

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

// TODO (fractalsbyx): This class seems to parallel DofManager quite a bit. Consider
// merging them?
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
                    const DofManager<dim>                      &_dof_manager,
                    const PDEOperatorBase<dim, degree, number> *_pde_operator);

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
  [[nodiscard]] std::vector<const dealii::AffineConstraints<number> *>
  get_change_constraints(const std::set<unsigned int> &field_indices,
                         unsigned int                  relative_level = 0) const;

  /**
   * @brief Getter function for the constraint of an index (constant reference).
   */
  [[nodiscard]] const dealii::AffineConstraints<number> &
  get_change_constraint(Types::Index index, unsigned int relative_level = 0) const;

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
  make_constraints_for_single_field(
    dealii::AffineConstraints<number>               &constraint,
    const dealii::DoFHandler<dim>                   &dof_handler,
    const std::map<unsigned int, BoundaryCondition> &bc_set,
    FieldInfo::TensorRank                            tensor_rank,
    Types::Index                                     field_index,
    bool                                             for_change_term);

  /**
   * @brief Add boundary conditions to a single constraint.
   */
  void
  make_bc_constraints(dealii::AffineConstraints<number>               &constraint,
                      const dealii::DoFHandler<dim>                   &dof_handler,
                      const std::map<unsigned int, BoundaryCondition> &boundary_condition,
                      FieldInfo::TensorRank                            tensor_rank,
                      Types::Index                                     field_index,
                      bool for_change_term = false);

  /**
   * @brief Apply constraints for common boundary conditions.
   */
  void
  make_one_boundary_constraint(dealii::AffineConstraints<number> &_constraints,
                               unsigned int                       boundary_id,
                               unsigned int                       component,
                               BoundaryCondition::Type            boundary_type,
                               number                             dirichlet_value,
                               const dealii::DoFHandler<dim>     &dof_handler,
                               FieldInfo::TensorRank              tensor_rank,
                               Types::Index                       field_index,
                               bool for_change_term = false) const;

  /**
   * @brief Apply natural constraints.
   */
  void
  make_natural_constraints() const;

  /**
   * @brief make dirichlet constraints.
   */
  void
  make_uniform_dirichlet_constraints(dealii::AffineConstraints<number> &_constraints,
                                     const dealii::DoFHandler<dim>     &dof_handler,
                                     const unsigned int                &boundary_id,
                                     const bool                        &is_vector_field,
                                     const number                      &value,

                                     const dealii::ComponentMask &mask) const;

  /**
   * @brief make nonuniform dirichlet constraints.
   */
  void
  make_nonuniform_dirichlet_constraints(dealii::AffineConstraints<number> &_constraints,
                                        const dealii::DoFHandler<dim>     &dof_handler,
                                        const unsigned int                &boundary_id,
                                        const unsigned int                &field_index,
                                        const bool                  &is_vector_field,
                                        const dealii::ComponentMask &mask,
                                        bool is_change_term = false) const;

  /**
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
                   FieldInfo::TensorRank              tensor_rank,
                   bool                               is_change_term = false) const;

  /**
   * @brief User-inputs.
   */
  const UserInputParameters<dim> *user_inputs;

  /**
   * @brief Dof manager pointer.
   */
  std::shared_ptr<const DofManager<dim>> dof_manager;

  /**
   * @brief PDE operator number.
   */
  const PDEOperatorBase<dim, degree, number> *pde_operator;

  /**
   * @brief Whether we have multigrid.
   */
  bool has_multigrid = false;

  /**
   * @brief Global minimum level for multigrid.
   */
  unsigned int global_min_level = 0;

  /**
   * @brief Constraints. Outer vector is indexed by field index. Inner vector is indexed
   * by relative mg level.
   */
  std::vector<std::vector<dealii::AffineConstraints<number>>> constraints;

  /**
   * @brief Constraints for Newton-Change solutions. Outer vector is indexed by field
   * index. Inner vector is indexed by relative mg level.
   */
  std::vector<std::vector<dealii::AffineConstraints<number>>> change_constraints;
};

PRISMS_PF_END_NAMESPACE
