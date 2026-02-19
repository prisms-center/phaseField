// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/field_container.h>
#include <prismspf/core/phase_field_tools.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class contains the user implementation of each PDE operator.
 */
template <unsigned int dim, unsigned int degree, typename number>
class PDEOperatorBase
{
public:
  using SizeType = dealii::VectorizedArray<number>;

  /**
   * @brief Constructor.
   */
  explicit PDEOperatorBase(const UserInputParameters<dim> &_user_inputs,
                           const PhaseFieldTools<dim>     &_pf_tools);

  /**
   * @brief Destructor.
   */
  virtual ~PDEOperatorBase() = default;

  /**
   * @brief User-implemented class for the setting initial conditions.
   */
  virtual void
  set_initial_condition(const unsigned int       &index,
                        const unsigned int       &component,
                        const dealii::Point<dim> &point,
                        number                   &scalar_value,
                        number                   &vector_component_value) const = 0;

  /**
   * @brief User-implemented class for the setting nonuniform boundary conditions.
   */
  virtual void
  set_nonuniform_dirichlet(const unsigned int       &index,
                           const unsigned int       &boundary_id,
                           const unsigned int       &component,
                           const dealii::Point<dim> &point,
                           number                   &scalar_value,
                           number                   &vector_component_value) const = 0;

  /**
   * @brief User-implemented class for the RHS of explicit equations.
   */
  virtual void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              unsigned int                         solver_id) const = 0;

  /**
   * @brief User-implemented class for the RHS of nonexplicit equations.
   */
  virtual void
  compute_lhs(FieldContainer<dim, degree, number> &variable_list,
              unsigned int                         solver_id) const = 0;

  /**
   * @brief Get the user inputs (constant reference).
   */
  [[nodiscard]] const UserInputParameters<dim> &
  get_user_inputs() const;

  /**
   * @brief Get the phase field tools (constant reference).
   */
  [[nodiscard]] const PhaseFieldTools<dim> &
  get_pf_tools() const;

private:
  /**
   * @brief The user-inputs.
   */
  const UserInputParameters<dim> *user_inputs;

  /**
   * @brief Phase field tools.
   */
  const PhaseFieldTools<dim> *pf_tools;
};

PRISMS_PF_END_NAMESPACE