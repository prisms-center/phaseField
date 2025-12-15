// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/affine_constraints.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/boundary_parameters.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

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
  ConstraintManager(
    const std::vector<FieldAttributes>                            &field_attributes,
    const std::set<SolveGroup>                                    &solve_groups,
    const UserInputParameters<dim>                                &_user_inputs,
    const std::shared_ptr<const PDEOperator<dim, degree, number>> &_pde_operator);

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
   * @brief Make constraints based on the inputs of the constructor.
   */
  void
  make_constraints(const dealii::Mapping<dim>                         &mapping,
                   const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers);

  /**
   * @brief Update time-dependent constraints.
   * For now this only updates the non-uniform dirichlet constraints.
   */
  void
  update_time_dependent_constraints(
    const dealii::Mapping<dim>                         &mapping,
    const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers);

private:
  /**
   * @brief Create a component mask.
   */
  [[nodiscard]] dealii::ComponentMask
  create_component_mask(unsigned int component, bool is_vector_field) const;

  /**
   * @brief Apply natural constraints.
   */
  void
  apply_natural_constraints() const;

  /**
   * @brief Apply dirichlet constraints.
   */
  template <typename num>
  void
  apply_dirichlet_constraints(const dealii::Mapping<dim>     &mapping,
                              const dealii::DoFHandler<dim>  &dof_handler,
                              const unsigned int             &boundary_id,
                              const bool                     &is_vector_field,
                              const num                      &value,
                              dealii::AffineConstraints<num> &_constraints,
                              const dealii::ComponentMask    &mask) const;

  /**
   * @brief Apply periodic constraints.
   */
  template <typename num>
  void
  apply_periodic_constraints(const dealii::DoFHandler<dim>  &dof_handler,
                             const unsigned int             &boundary_id,
                             dealii::AffineConstraints<num> &_constraints,
                             const dealii::ComponentMask    &mask) const;

  /**
   * @brief Apply nonuniform dirichlet constraints.
   */
  template <typename num>
  void
  apply_nonuniform_dirichlet_constraints(const dealii::Mapping<dim>     &mapping,
                                         const dealii::DoFHandler<dim>  &dof_handler,
                                         const unsigned int             &boundary_id,
                                         const unsigned int             &index,
                                         const bool                     &is_vector_field,
                                         dealii::AffineConstraints<num> &_constraints,
                                         const dealii::ComponentMask    &mask,
                                         bool is_change_term = false) const;

  /**
   * @brief Make the constraint for a single index.
   */
  void
  make_constraint(const dealii::Mapping<dim>    &mapping,
                  const dealii::DoFHandler<dim> &dof_handler,
                  unsigned int                   index);

  /**
   * @brief Set the dirichlet constraint for the pinned point.
   */
  template <typename num>
  void
  set_pinned_point(const dealii::DoFHandler<dim>  &dof_handler,
                   dealii::AffineConstraints<num> &_constraints,
                   unsigned int                    index,
                   bool                            is_change_term = false) const;

  /**
   * @brief Clear, reinitialize and make hanging node constraints
   */
  template <typename num>
  void
  apply_generic_constraints(const dealii::DoFHandler<dim>  &dof_handler,
                            dealii::AffineConstraints<num> &_constraints) const;

  /**
   * @brief Apply constraints for common boundary conditions.
   */
  template <typename num, int spacedim>
  void
  apply_constraints(const dealii::Mapping<dim>     &mapping,
                    const dealii::DoFHandler<dim>  &dof_handler,
                    dealii::AffineConstraints<num> &_constraints,
                    const BoundaryCondition        &boundary_condition,
                    BoundaryCondition::Type         boundary_type,
                    unsigned int                    boundary_id,
                    unsigned int                    component,
                    unsigned int                    index,
                    bool                            is_change_term = false) const;

  /**
   * @brief User-inputs.
   */
  const UserInputParameters<dim> *user_inputs;

  /**
   * @brief PDE operator number.
   */
  std::shared_ptr<const PDEOperator<dim, degree, number>> pde_operator;

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
  std::vector<std::vector<std::shared_ptr<dealii::AffineConstraints<float>>>> constraints;
};

PRISMS_PF_END_NAMESPACE
