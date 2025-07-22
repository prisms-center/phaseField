// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/concurrent_solver.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class handles the explicit solves of all constant fields
 */
template <unsigned int dim, unsigned int degree, typename number>
class ConcurrentConstantSolver : public ConcurrentSolver<dim, degree, number>
{
public:
  /**
   * @brief Constructor.
   */
  explicit ConcurrentConstantSolver(const SolverContext<dim, degree> &_solver_context,
                                    Types::Index                      _solve_priority = 0)
    : ConcurrentSolver<dim, degree, number>(_solver_context,
                                            FieldSolveType::ExplicitConstant,
                                            _solve_priority) {};

  /**
   * @brief Destructor.
   */
  ~ConcurrentConstantSolver() override = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  ConcurrentConstantSolver(const ConcurrentConstantSolver &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  ConcurrentConstantSolver &
  operator=(const ConcurrentConstantSolver &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  ConcurrentConstantSolver(ConcurrentConstantSolver &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  ConcurrentConstantSolver &
  operator=(ConcurrentConstantSolver &&solver_base) noexcept = delete;

  /**
   * @brief Initialize the solver.
   */
  void
  init() override
  {
    // Call the base class init
    this->ConcurrentSolver<dim, degree, number>::init();

    // Do nothing
  };

  /**
   * @brief Reinitialize the solver.
   */
  void
  reinit() override
  {
    // Call the base class reinit
    this->ConcurrentSolver<dim, degree, number>::reinit();

    // Do nothing
  };

  /**
   * @brief Solve for a single update step.
   */
  void
  solve() override
  {
    // Call the base class solve
    this->ConcurrentSolver<dim, degree, number>::solve();

    // Do nothing
  };

  /**
   * @brief Print information about the solver to summary.log.
   */
  void
  print()
  {
    // Print the base class information
    this->ConcurrentSolver<dim, degree, number>::print();
  }
};

PRISMS_PF_END_NAMESPACE