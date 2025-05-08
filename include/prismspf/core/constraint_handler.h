// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/affine_constraints.h>

#include <prismspf/core/multigrid_info.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class userInputParameters;

/**
 * \brief The class handles the generation and application of boundary conditions based on
 * the user-inputs.
 */
template <unsigned int dim>
class constraintHandler
{
public:
  /**
   * \brief Constructor.
   */
  constraintHandler(const userInputParameters<dim> &_user_inputs,
                    const MGInfo<dim>              &_mg_info);

  /**
   * \brief Getter function for the constraints.
   */
  [[nodiscard]] std::vector<const dealii::AffineConstraints<double> *>
  get_constraints();

  /**
   * \brief Getter function for the constraint of an index (constant reference).
   */
  [[nodiscard]] const dealii::AffineConstraints<double> &
  get_constraint(unsigned int index) const;

  /**
   * \brief Getter function for the multigrid constraints of a certain level.
   */
  [[nodiscard]] std::vector<const dealii::AffineConstraints<float> *>
  get_mg_constraints(unsigned int level);

  /**
   * \brief Getter function for the multigrid constraint of a certain level and index
   * (constant reference).
   */
  [[nodiscard]] const dealii::AffineConstraints<float> &
  get_mg_constraint(unsigned int level, unsigned int index) const;

  /**
   * \brief Make constraints based on the inputs of the constructor.
   */
  void
  make_constraints(const dealii::Mapping<dim>                         &mapping,
                   const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers);

  /**
   * \brief Make multigrid constraints for a given level based on the inputs of the
   * constructor.
   */
  void
  make_mg_constraints(const dealii::Mapping<dim>                         &mapping,
                      const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers,
                      unsigned int                                        level);

  /**
   * \brief Update time-dependent constraints.
   *
   * For now this only updates the non-uniform dirichlet constraints.
   */
  void
  update_time_dependent_constraints(
    const dealii::Mapping<dim>                         &mapping,
    const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers);

  /**
   * \brief Update time-dependent multigrid constraints for a given level.
   *
   * For now this only updates the non-uniform dirichlet constraints.
   */
  void
  update_time_dependent_mg_constraints(
    const dealii::Mapping<dim>                         &mapping,
    const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers,
    unsigned int                                        level);

private:
  /**
   * \brief Make the constraint for a single index.
   */
  void
  make_constraint(const dealii::Mapping<dim>    &mapping,
                  const dealii::DoFHandler<dim> &dof_handler,
                  unsigned int                   index);

  /**
   * \brief Make the multigrid constraint for a single index at a single level.
   */
  void
  make_mg_constraint(const dealii::Mapping<dim>    &mapping,
                     const dealii::DoFHandler<dim> &dof_handler,
                     unsigned int                   index,
                     unsigned int                   level,
                     dependencyType                 dependency_type);

  /**
   * \brief Set the dirichlet constraint for the pinned point.
   */
  template <typename number>
  void
  set_pinned_point(const dealii::DoFHandler<dim>     &dof_handler,
                   dealii::AffineConstraints<number> &constraints,
                   unsigned int                       index) const;

  /**
   * \brief Set the dirichlet constraint for the pinned point.
   */
  template <typename number>
  void
  set_mg_pinned_point(const dealii::DoFHandler<dim>     &dof_handler,
                      dealii::AffineConstraints<number> &constraints,
                      unsigned int                       index) const;

  /**
   * \brief Clear, reinitialize and make hanging node constraints
   */
  template <typename number>
  void
  apply_generic_constraints(const dealii::DoFHandler<dim>     &dof_handler,
                            dealii::AffineConstraints<number> &constraints) const;

  /**
   * \brief Apply constraints for common boundary conditions.
   */
  template <typename number, int spacedim>
  void
  apply_constraints(const dealii::Mapping<dim>        &mapping,
                    const dealii::DoFHandler<dim>     &dof_handler,
                    dealii::AffineConstraints<number> &constraints,
                    const boundaryCondition           &boundary_condition,
                    boundaryCondition::type            boundary_type,
                    unsigned int                       boundary_id,
                    unsigned int                       component,
                    unsigned int                       index) const;

  /**
   * \brief Apply multigrid constraints for common boundary conditions. The only
   * difference between this function and the previous one is that an dirichlet boundary
   * conditions are constrained to 0.
   */
  template <typename number, int spacedim>
  void
  apply_mg_constraints(const dealii::Mapping<dim>        &mapping,
                       const dealii::DoFHandler<dim>     &dof_handler,
                       dealii::AffineConstraints<number> &constraints,
                       boundaryCondition::type            boundary_type,
                       unsigned int                       boundary_id,
                       unsigned int                       component) const;

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> *user_inputs;

  /**
   * \brief Multigrid info
   */
  const MGInfo<dim> *mg_info;

  /**
   * \brief Whether we have multigrid.
   */
  bool has_multigrid = false;

  /**
   * \brief Global minimum level for multigrid.
   */
  unsigned int global_min_level = 0;

  /**
   * \brief Constraints.
   */
  std::vector<dealii::AffineConstraints<double>> constraints;

  /**
   * \brief Multigrid constraints.
   */
  std::vector<std::vector<dealii::AffineConstraints<float>>> mg_constraints;
};

PRISMS_PF_END_NAMESPACE
