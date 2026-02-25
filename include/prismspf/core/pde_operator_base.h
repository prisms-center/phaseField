// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/field_container.h>
#include <prismspf/core/phase_field_tools.h>
#include <prismspf/core/simulation_timer.h>
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
                           const PhaseFieldTools<dim>     &_pf_tools)
    : user_inputs(&_user_inputs)
    , pf_tools(&_pf_tools)
  {}

  /**
   * @brief Destructor.
   */
  virtual ~PDEOperatorBase() = default;

  /**
   * @brief User-implemented class for the setting initial conditions.
   */
  virtual void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const
  {}

  /**
   * @brief User-implemented class for the setting nonuniform boundary conditions.
   */
  virtual void
  set_nonuniform_dirichlet([[maybe_unused]] const unsigned int       &index,
                           [[maybe_unused]] const unsigned int       &boundary_id,
                           [[maybe_unused]] const unsigned int       &component,
                           [[maybe_unused]] const dealii::Point<dim> &point,
                           [[maybe_unused]] number                   &scalar_value,
                           [[maybe_unused]] number &vector_component_value) const
  {}

  /**
   * @brief User-implemented class for the RHS of explicit equations.
   */
  virtual void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int                         solver_id) const
  {}

  /**
   * @brief User-implemented class for the RHS of nonexplicit equations.
   */
  virtual void
  compute_lhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int                         solver_id) const
  {}

  /**
   * @brief Get the user inputs (constant reference).
   */
  [[nodiscard]] const UserInputParameters<dim> &
  get_user_inputs() const
  {
    return *user_inputs;
  }

  /**
   * @brief Get the phase field tools (constant reference).
   */
  [[nodiscard]] const PhaseFieldTools<dim> &
  get_pf_tools() const
  {
    return *pf_tools;
  }

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