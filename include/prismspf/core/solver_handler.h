// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/types.h>

#include <prismspf/solvers/concurrent_constant_solver.h>
#include <prismspf/solvers/concurrent_explicit_postprocess_solver.h>
#include <prismspf/solvers/concurrent_explicit_solver.h>
#include <prismspf/solvers/sequential_auxiliary_solver.h>
#include <prismspf/solvers/sequential_co_nonlinear_solver.h>
#include <prismspf/solvers/sequential_linear_solver.h>
#include <prismspf/solvers/sequential_self_nonlinear_solver.h>

#include <prismspf/config.h>

#include <map>
#include <set>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolverContext;

template <unsigned int dim, unsigned int degree, typename number>
class SolveBlock
{
public:
  /**
   * @brief Constructor.
   */
  explicit SolveBlock(const SolverContext<dim, degree, number> &_solver_context,
                      Types::Index                              block_index);

  /**
   * @brief Destructor.
   */
  ~SolveBlock() = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  SolveBlock(const SolveBlock &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  SolveBlock &
  operator=(const SolveBlock &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  SolveBlock(SolveBlock &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  SolveBlock &
  operator=(SolveBlock &&solver_base) noexcept = delete;

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
   * @brief Explicit constant field solver class.
   */
  ConcurrentConstantSolver<dim, degree, number> concurrent_constant_solver;

  /**
   * @brief Explicit field solver class.
   */
  ConcurrentExplicitSolver<dim, degree, number> concurrent_explicit_solver;

  /**
   * @brief Postprocessed explicit field solver class.
   */
  ConcurrentExplicitPostprocessSolver<dim, degree, number>
    concurrent_explicit_postprocess_solver;

  /**
   * @brief Nonexplicit auxiliary field solver class.
   */
  SequentialAuxiliarySolver<dim, degree, number> sequential_auxiliary_solver;

  /**
   * @brief Nonexplicit linear field solver class.
   */
  SequentialLinearSolver<dim, degree, number> sequential_linear_solver;

  /**
   * @brief Nonexplicit self-nonlinear field solver class.
   */
  SequentialSelfNonlinearSolver<dim, degree, number> sequential_self_nonlinear_solver;

  /**
   * @brief Nonexplicit co-nonlinear field solver class.
   */
  SequentialCoNonlinearSolver<dim, degree, number> sequential_co_nonlinear_solver;
};

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
  std::map<Types::Index, SolveBlock<dim, degree, number>> solve_blocks;
};

PRISMS_PF_END_NAMESPACE