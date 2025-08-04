// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/solver_base.h>

#include <prismspf/config.h>

#include <functional>
#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolverContext;

template <unsigned int dim, unsigned int degree, typename number>
class ConcurrentSolver : public SolverBase<dim, degree, number>
{
public:
  /**
   * @brief Constructor.
   */
  ConcurrentSolver(const SolverContext<dim, degree, number> &_solver_context,
                   const FieldSolveType                     &_field_solve_type,
                   Types::Index                              _solve_priority = 0);

  /**
   * @brief Destructor.
   */
  ~ConcurrentSolver() override = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  ConcurrentSolver(const ConcurrentSolver &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  ConcurrentSolver &
  operator=(const ConcurrentSolver &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  ConcurrentSolver(ConcurrentSolver &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  ConcurrentSolver &
  operator=(ConcurrentSolver &&solver_base) noexcept = delete;

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

  /**
   * @brief Solve the explicit equations
   *
   * This is a common function for solving explicit (RHS only) equations that are
   * independent of one another.
   *
   * Rather than duplicate this code a bunch of times for explicit, postprocess, amr,
   * nucleation, etc... fields we have it here. Importantly, some of these have different
   * functions, so we require this function.
   */
  void
  solve_explicit_equations(
    const std::function<
      void(std::vector<typename SolverBase<dim, degree, number>::VectorType *> &,
           const std::vector<typename SolverBase<dim, degree, number>::VectorType *> &)>
      &function);

  /**
   * @brief Get the system matrix.
   */
  [[nodiscard]] std::unique_ptr<
    typename SolverBase<dim, degree, number>::SystemMatrixType> &
  get_system_matrix()
  {
    return system_matrix;
  }

  /**
   * @brief Get the mapping from global solution vectors to the local ones.
   */
  [[nodiscard]] const std::vector<Types::Index> &
  get_global_to_local_solution_mapping()
  {
    return global_to_local_solution;
  }

  /**
   * @brief Get the src solution subset.
   */
  [[nodiscard]] const std::vector<
    typename SolverBase<dim, degree, number>::VectorType *> &
  get_src_solution_subset()
  {
    return solution_subset;
  }

  /**
   * @brief Get the dst solution subset.
   */
  [[nodiscard]] std::vector<typename SolverBase<dim, degree, number>::VectorType *> &
  get_dst_solution_subset()
  {
    return new_solution_subset;
  }

private:
  /**
   * @brief Matrix-free operator.
   */
  std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>
    system_matrix;

  /**
   * @brief Mapping from global solution vectors to the local ones
   */
  std::vector<Types::Index> global_to_local_solution;

  /**
   * @brief Subset of solutions fields that are necessary for concurrent solves.
   */
  std::vector<typename SolverBase<dim, degree, number>::VectorType *> solution_subset;

  /**
   * @brief Subset of new solutions fields that are necessary for concurrent solves.
   */
  std::vector<typename SolverBase<dim, degree, number>::VectorType *> new_solution_subset;
};

PRISMS_PF_END_NAMESPACE
