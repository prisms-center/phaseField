// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/point.h>

#include <prismspf/core/types.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

template <unsigned int dim, unsigned int degree, typename number>
class VariableContainer;

template <unsigned int dim>
struct PhaseFieldTools;

/**
 * @brief This class contains the user implementation of each PDE operator.
 */
template <unsigned int dim, unsigned int degree, typename number>
class PDEOperator
{
public:
  using SizeType = dealii::VectorizedArray<number>;

  /**
   * @brief Constructor.
   */
  explicit PDEOperator(const UserInputParameters<dim> &_user_inputs,
                       const PhaseFieldTools<dim>     &_pf_tools);

  /**
   * @brief Destructor.
   */
  virtual ~PDEOperator() = default;

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
  compute_explicit_rhs(VariableContainer<dim, degree, number> &variable_list,
                       const dealii::Point<dim, SizeType>     &q_point_loc,
                       const SizeType                         &element_volume,
                       Types::Index                            solve_block) const = 0;

  /**
   * @brief User-implemented class for the RHS of nonexplicit equations.
   */
  virtual void
  compute_nonexplicit_rhs(VariableContainer<dim, degree, number> &variable_list,
                          const dealii::Point<dim, SizeType>     &q_point_loc,
                          const SizeType                         &element_volume,
                          Types::Index                            solve_block,
                          Types::Index index = Numbers::invalid_index) const = 0;

  /**
   * @brief User-implemented class for the LHS of nonexplicit equations.
   */
  virtual void
  compute_nonexplicit_lhs(VariableContainer<dim, degree, number> &variable_list,
                          const dealii::Point<dim, SizeType>     &q_point_loc,
                          const SizeType                         &element_volume,
                          Types::Index                            solve_block,
                          Types::Index index = Numbers::invalid_index) const = 0;

  /**
   * @brief User-implemented class for the RHS of postprocessed explicit equations.
   */
  virtual void
  compute_postprocess_explicit_rhs(VariableContainer<dim, degree, number> &variable_list,
                                   const dealii::Point<dim, SizeType>     &q_point_loc,
                                   const SizeType                         &element_volume,
                                   Types::Index solve_block) const = 0;

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

  /**
   * @brief Get the timestep (copy).
   */
  [[nodiscard]] number
  get_timestep() const;

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