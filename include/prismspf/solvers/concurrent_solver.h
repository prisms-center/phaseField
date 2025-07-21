// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/solver_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class ConcurrentSolver : public SolverBase<dim, degree, number>
{
public:
  /**
   * @brief Constructor.
   */
  ConcurrentSolver(const SolverContext<dim, degree> &_solver_context,
                   const FieldSolveType             &_field_solve_type,
                   unsigned int                      _solve_priority = 0);

  /**
   * @brief Destructor.
   */
  virtual ~ConcurrentSolver() = 0;

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
  init() override
  {
    // Call the base class init
    this->SolverBase<dim, degree, number>::init();

    // Compute the shared dependencies
    this->compute_shared_dependencies(this->field_solve_type, this->solve_priority);
  };

  /**
   * @brief Reinitialize the solver.
   */
  void
  reinit() override
  {
    // Call the base class reinit
    this->SolverBase<dim, degree, number>::reinit();
  };

  /**
   * @brief Solve for a single update step.
   */
  void
  solve() override
  {
    // Call the base class solve
    this->SolverBase<dim, degree, number>::solve();
  };

  /**
   * @brief Print information about the solver to summary.log.
   */
  void
  print()
  {
    // Print the base class information
    this->SolverBase<dim, degree, number>::print();
  }

  /**
   * @brief Get the system matrix.
   */
  [[nodiscard]] std::unique_ptr<
    typename SolverBase<dim, degree, number>::SystemMatrixType> &
  get_system_matrix()
  {
    return system_matrix;
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
  std::vector<std::vector<Types::Index>> global_to_local_solution;

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
