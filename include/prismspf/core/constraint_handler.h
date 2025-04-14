// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/affine_constraints.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
class userInputParameters;

/**
 * \brief The class handles the generation and application of boundary conditions based on
 * the user-inputs.
 */
template <int dim>
class constraintHandler
{
public:
  /**
   * \brief Constructor.
   */
  explicit constraintHandler(const userInputParameters<dim> &_user_inputs);

  /**
   * \brief Getter function for the constraints.
   */
  [[nodiscard]] std::vector<const dealii::AffineConstraints<double> *>
  get_constraints();

  /**
   * \brief Getter function for the constraint of an index (constant reference).
   */
  [[nodiscard]] const dealii::AffineConstraints<double> &
  get_constraint(const unsigned int &index) const;

  /**
   * \brief Getter function for the multigrid constraints of an index (constant
   * reference).
   */
  [[nodiscard]] const dealii::MGLevelObject<dealii::AffineConstraints<float>> &
  get_mg_constraint(unsigned int index) const;

  /**
   * \brief Getter function for the multigrid constraints of an index at a certain
   * multigrid level (constant reference).
   */
  [[nodiscard]] const dealii::AffineConstraints<float> &
  get_mg_constraint(unsigned int index, unsigned int level) const;

  /**
   * \brief Getter function for the multigrid constraints at a certain multigrid level
   * (constant reference).
   */
  [[nodiscard]] std::vector<dealii::AffineConstraints<float> *>
  get_mg_level_constraints(unsigned int level) const;

  /**
   * \brief Make constraints based on the inputs of the constructor.
   */
  void
  make_constraints(const dealii::Mapping<dim>                         &mapping,
                   const std::vector<const dealii::DoFHandler<dim> *> &dof_handlers);

  /**
   * \brief Make multigrid constraints based on the inputs of the constructor.
   */
  void
  make_mg_constraints(
    const dealii::Mapping<dim> &mapping,
    const std::map<unsigned int, dealii::MGLevelObject<dealii::DoFHandler<dim>>>
      &mg_dof_handlers);

private:
  /**
   * \brief Make the constrainst for a single index.
   */
  void
  make_constraint(const dealii::Mapping<dim>    &mapping,
                  const dealii::DoFHandler<dim> &dof_handler,
                  const unsigned int            &index);

  /**
   * \brief Make the multigrid constrainst for a single index.
   */
  void
  make_mg_constraint(const dealii::Mapping<dim>                           &mapping,
                     const dealii::MGLevelObject<dealii::DoFHandler<dim>> &dof_handler,
                     const unsigned int                                   &index);

  /**
   * \brief Set the dirichlet constraint for the pinned point.
   */
  void
  set_pinned_point(const dealii::DoFHandler<dim> &dof_handler, const unsigned int &index);

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> *user_inputs;

  /**
   * \brief Constraints.
   */
  std::vector<dealii::AffineConstraints<double>> constraints;

  /**
   * \brief Multigrid constraints.
   */
  std::map<unsigned int, dealii::MGLevelObject<dealii::AffineConstraints<float>>>
    mg_constraints;
};

PRISMS_PF_END_NAMESPACE
