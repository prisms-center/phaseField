// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/types.h>

#include <prismspf/solvers/concurrent_solver.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolverContext;

/**
 * @brief This class handles the explicit solves of all explicit fields
 */
template <unsigned int dim, unsigned int degree, typename number>
class ConcurrentExplicitSolver : public ConcurrentSolver<dim, degree, number>
{
public:
  /**
   * @brief Constructor.
   */
  explicit ConcurrentExplicitSolver(
    const SolverContext<dim, degree, number> &_solver_context,
    Types::Index                              _solve_priority = 0);

  /**
   * @brief Destructor.
   */
  ~ConcurrentExplicitSolver() override = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  ConcurrentExplicitSolver(const ConcurrentExplicitSolver &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  ConcurrentExplicitSolver &
  operator=(const ConcurrentExplicitSolver &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  ConcurrentExplicitSolver(ConcurrentExplicitSolver &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  ConcurrentExplicitSolver &
  operator=(ConcurrentExplicitSolver &&solver_base) noexcept = delete;

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
  solve() override;

  /**
   * @brief Print information about the solver to summary.log.
   */
  void
  print() override;
};

PRISMS_PF_END_NAMESPACE
