// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/sequential_solver.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class handles all linear solves.
 */
template <unsigned int dim, unsigned int degree, typename number>
class SequentialLinearSolver : public SequentialSolver<dim, degree, number>
{
public:
  /**
   * @brief Constructor.
   */
  explicit SequentialLinearSolver(const SolverContext<dim, degree> &_solver_context,
                                  Types::Index                      _solve_priority = 0)
    : SequentialSolver<dim, degree, number>(_solver_context,
                                            FieldSolveType::NonexplicitLinear,
                                            _solve_priority) {};

  /**
   * @brief Destructor.
   */
  ~SequentialLinearSolver() override = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialLinearSolver(const SequentialLinearSolver &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialLinearSolver &
  operator=(const SequentialLinearSolver &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialLinearSolver(SequentialLinearSolver &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialLinearSolver &
  operator=(SequentialLinearSolver &&solver_base) noexcept = delete;

  /**
   * @brief Initialize the solver.
   */
  void
  init() override
  {
    // Call the base class init
    this->SequentialSolver<dim, degree, number>::init();

    // If the solver is empty we can just return early.
    if (this->solver_is_empty())
      {
        return;
      }

    // Init the linear solvers
    for (const auto &[index, variable] : this->get_subset_attributes())
      {
        this->init_linear_solver(variable);
      }
  };

  /**
   * @brief Reinitialize the solver.
   */
  void
  reinit() override
  {
    // Call the base class reinit
    this->SequentialSolver<dim, degree, number>::reinit();

    // If the solver is empty we can just return early.
    if (this->solver_is_empty())
      {
        return;
      }

    // Reinit the linear solvers
    for (const auto &[index, variable] : this->get_subset_attributes())
      {
        this->reinit_linear_solver(variable);
      }
  };

  /**
   * @brief Solve for a single update step.
   */
  void
  solve() override
  {
    // Call the base class solve
    this->SequentialSolver<dim, degree, number>::solve();

    // If the solver is empty we can just return early.
    if (this->solver_is_empty())
      {
        return;
      }

    // Solve each field
    for (const auto &[index, variable] : this->get_subset_attributes())
      {
        this->solve_linear_solver(variable);
      }
  };

  /**
   * @brief Print information about the solver to summary.log.
   */
  void
  print()
  {
    // Print the base class information
    this->SequentialSolver<dim, degree, number>::print();
  }
};

PRISMS_PF_END_NAMESPACE
