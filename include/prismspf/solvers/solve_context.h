// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/constraint_manager.h>
#include <prismspf/core/dof_manager.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_indexer.h>
#include <prismspf/core/triangulation_manager.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class provides context for a solver with ptrs to all the relevant
 * dependencies.
 *
 * The context in this case refers to all the finite element machinery needed to solve the
 * fields. For example, it contains the triangulation manager and pde operator that
 * evaluates the user-specified PDEs.
 */
template <unsigned int dim, unsigned int degree, typename number>
class SolveContext
{
public:
  /**
   * @brief Constructor.
   */
  SolveContext(const UserInputParameters<dim>               &_user_inputs,
               const TriangulationManager<dim>              &_triangulation_manager,
               const ConstraintManager<dim, degree, number> &_constraint_manager,
               const DofManager<dim>                        &_dof_manager,
               SolutionIndexer<dim, number>                 &_solution_indexer,
               std::shared_ptr<const PDEOperator<dim, degree, number>> _pde_operator)
    : user_inputs(&_user_inputs)
    , triangulation_manager(&_triangulation_manager)
    , constraint_manager(&_constraint_manager)
    , dof_manager(&_dof_manager)
    , solution_indexer(&_solution_indexer)
    , pde_operator(_pde_operator) {};

  /**
   * @brief Get the user-inputs.
   */
  [[nodiscard]] const UserInputParameters<dim> &
  get_user_inputs() const
  {
    Assert(user_inputs != nullptr, dealii::ExcNotInitialized());
    return *user_inputs;
  }

  /**
   * @brief Get the triangulation manager.
   */
  [[nodiscard]] const TriangulationManager<dim> &
  get_triangulation_manager() const
  {
    Assert(triangulation_manager != nullptr, dealii::ExcNotInitialized());
    return *triangulation_manager;
  }

  /**
   * @brief Get the constraint manager.
   */
  [[nodiscard]] const ConstraintManager<dim, degree, number> &
  get_constraint_manager() const
  {
    Assert(constraint_manager != nullptr, dealii::ExcNotInitialized());
    return *constraint_manager;
  }

  /**
   * @brief Get the dof manager.
   */
  [[nodiscard]] const DofManager<dim> &
  get_dof_manager() const
  {
    Assert(dof_manager != nullptr, dealii::ExcNotInitialized());
    return *dof_manager;
  }

  /**
   * @brief Get the solution manager.
   */
  [[nodiscard]] SolutionIndexer<dim, number> &
  get_solution_indexer() const
  {
    Assert(solution_indexer != nullptr, dealii::ExcNotInitialized());
    return *solution_indexer;
  }

  /**
   * @brief Get a shared pointer to the pde operator.
   */
  [[nodiscard]] const std::shared_ptr<const PDEOperator<dim, degree, number>> &
  get_pde_operator() const
  {
    Assert(pde_operator != nullptr, dealii::ExcNotInitialized());
    return pde_operator;
  }

private:
  /**
   * @brief User-inputs.
   */
  const UserInputParameters<dim> *user_inputs;

  /**
   * @brief Triangulation manager.
   */
  const TriangulationManager<dim> *triangulation_manager;

  /**
   * @brief Constraint manager.
   */
  const ConstraintManager<dim, degree, number> *constraint_manager;

  /**
   * @brief DoF manager.
   */
  const DofManager<dim> *dof_manager;

  /**
   * @brief Solution manager.
   */
  SolutionIndexer<dim, number> *solution_indexer;

  /**
   * @brief PDE operator.
   */
  std::shared_ptr<const PDEOperator<dim, degree, number>> pde_operator;

  /**
   * @brief PDE operator for float precision.
   */
  std::shared_ptr<const PDEOperator<dim, degree, float>> pde_operator_float;
};

PRISMS_PF_END_NAMESPACE
