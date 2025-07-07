// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/solver_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SequentialSolver : public SolverBase<dim, degree, number>
{
public:
  /**
   * \brief Constructor.
   */
  SequentialSolver(const SolverContext<dim, degree> &_solver_context,
                   const FieldSolveType             &_field_solve_type,
                   unsigned int                      _solve_priority = 0);

  /**
   * \brief Destructor.
   */
  virtual ~SequentialSolver() = 0;

  /**
   * \brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialSolver(const SequentialSolver &solver_base) = delete;

  /**
   * \brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialSolver &
  operator=(const SequentialSolver &solver_base) = delete;

  /**
   * \brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialSolver(SequentialSolver &&solver_base) noexcept = delete;

  /**
   * \brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialSolver &
  operator=(SequentialSolver &&solver_base) noexcept = delete;

  /**
   * \brief Initialize the solver.
   */
  void
  init() override
  {
    // Call the base class init
    this->SolverBase<dim, degree, number>::init();
  };

  /**
   * \brief Reinitialize the solver.
   */
  void
  reinit() override
  {
    // Call the base class reinit
    this->SolverBase<dim, degree, number>::reinit();
  };

  /**
   * \brief Solve for a single update step.
   */
  void
  solve() override
  {
    // Call the base class solve
    this->SolverBase<dim, degree, number>::solve();
  };

  /**
   * \brief Print information about the solver to summary.log.
   */
  void
  print()
  { // Print the base class information
    this->SolverBase<dim, degree, number>::print();
  };

  /**
   * \brief Get the matrix-free operator for the residual side.
   */
  [[nodiscard]] std::map<
    unsigned int,
    std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>> &
  get_system_matrix()
  {
    return system_matrix;
  }

  /**
   * \brief Get the matrix-free operator for the newton update side.
   */
  [[nodiscard]] std::map<
    unsigned int,
    std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>> &
  get_update_system_matrix()
  {
    return update_system_matrix;
  }

private:
  /**
   * \brief Matrix-free operator for the residual side.
   */
  std::map<unsigned int,
           std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>>
    system_matrix;

  /**
   * \brief Matrix-free operator for the newton update side.
   */
  std::map<unsigned int,
           std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>>
    update_system_matrix;
};

PRISMS_PF_END_NAMESPACE
