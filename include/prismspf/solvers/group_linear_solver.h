// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/types.h>

#include <prismspf/solvers/group_solver_base.h>
#include <prismspf/solvers/mf_operator.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolverContext;

/**
 * @brief This class handles the explicit solves of all explicit fields
 */
template <unsigned int dim, unsigned int degree, typename number>
class LinearSolver : public GroupSolverBase<dim, degree, number>
{
public:
  using GroupSolverBase = GroupSolverBase<dim, degree, number>;
  using GroupSolverBase::mf_operators;
  using GroupSolverBase::solutions;
  using GroupSolverBase::solve_group;
  using GroupSolverBase::solver_context;

  /**
   * @brief Constructor.
   */
  LinearSolver(SolveGroup                                _solve_group,
               const SolverContext<dim, degree, number> &_solver_context)
    : GroupSolverBase(_solve_group, _solver_context)
  {}

  /**
   * @brief Destructor.
   */
  ~LinearSolver() override = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  LinearSolver(const LinearSolver &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  LinearSolver &
  operator=(const LinearSolver &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  LinearSolver(LinearSolver &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  LinearSolver &
  operator=(LinearSolver &&solver_base) noexcept = delete;

  /**
   * @brief Initialize the solver.
   */
  void
  init() override;

  /**
   * @brief Reinitialize the solver.
   */
  void
  reinit() override;

  /**
   * @brief Solve for a single update step.
   */
  void
  solve_level(unsigned int relative_level) override;

  /**
   * @brief Solve for a single update step.
   */
  void
  solve() override;

  /**
   * @brief Print information about the solver to summary.log.
   */
  void
  print() override;

private:
  /**
   * @brief Matrix free operators for each level
   */
  std::vector<MFOperator<dim, degree, number>> lhs_mf_operators;
};

PRISMS_PF_END_NAMESPACE
