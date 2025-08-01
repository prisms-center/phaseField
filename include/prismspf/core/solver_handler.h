// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/concurrent_constant_solver.h>
#include <prismspf/solvers/concurrent_explicit_postprocess_solver.h>
#include <prismspf/solvers/concurrent_explicit_solver.h>
#include <prismspf/solvers/sequential_auxiliary_solver.h>
#include <prismspf/solvers/sequential_co_nonlinear_solver.h>
#include <prismspf/solvers/sequential_linear_solver.h>
#include <prismspf/solvers/sequential_self_nonlinear_solver.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief The class handles the initialization and solving of the various types of solvers
 * in PRISMS-PF.
 */
template <unsigned int dim, unsigned int degree, typename number>
class SolverHandler
{
public:
  /**
   * @brief Constructor.
   */
  explicit SolverHandler(const SolverContext<dim, degree, number> &_solver_context);

  /**
   * @brief Destructor.
   */
  ~SolverHandler() = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  SolverHandler(const SolverHandler &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  SolverHandler &
  operator=(const SolverHandler &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  SolverHandler(SolverHandler &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  SolverHandler &
  operator=(SolverHandler &&solver_base) noexcept = delete;

  /**
   * @brief Initialize all of the solvers.
   */
  void
  init();

  /**
   * @brief Reinitialize all of the solvers.
   */
  void
  reinit();

  /**
   * @brief Solve a single incremenet for all of the solvers.
   *
   * TODO (landinjm): Fix the logic here with the passed variables
   */
  void
  solve(unsigned int increment, bool update_postprocessed);

private:
  /**
   * @brief Set of solve blocks that we have.
   */
  std::set<Types::Index> solve_blocks;

  /**
   * @brief Explicit constant field solver class.
   */
  std::map<Types::Index, ConcurrentConstantSolver<dim, degree, number>>
    concurrent_constant_solver;

  /**
   * @brief Explicit field solver class.
   */
  std::map<Types::Index, ConcurrentExplicitSolver<dim, degree, number>>
    concurrent_explicit_solver;

  /**
   * @brief Postprocessed explicit field solver class.
   */
  std::map<Types::Index, ConcurrentExplicitPostprocessSolver<dim, degree, number>>
    concurrent_explicit_postprocess_solver;

  /**
   * @brief Nonexplicit auxiliary field solver class.
   */
  std::map<Types::Index, SequentialAuxiliarySolver<dim, degree, number>>
    sequential_auxiliary_solver;

  /**
   * @brief Nonexplicit linear field solver class.
   */
  std::map<Types::Index, SequentialLinearSolver<dim, degree, number>>
    sequential_linear_solver;

  /**
   * @brief Nonexplicit self-nonlinear field solver class.
   */
  std::map<Types::Index, SequentialSelfNonlinearSolver<dim, degree, number>>
    sequential_self_nonlinear_solver;

  /**
   * @brief Nonexplicit co-nonlinear field solver class.
   */
  std::map<Types::Index, SequentialCoNonlinearSolver<dim, degree, number>>
    sequential_co_nonlinear_solver;
};

PRISMS_PF_END_NAMESPACE