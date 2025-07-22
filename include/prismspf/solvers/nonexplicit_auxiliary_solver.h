// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/timer.h>

#include <prismspf/solvers/sequential_solver.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class handles all auxiliary solves.
 */
template <unsigned int dim, unsigned int degree, typename number>
class SequentialAuxiliarySolver : public SequentialSolver<dim, degree, number>
{
public:
  /**
   * @brief Constructor.
   */
  explicit SequentialAuxiliarySolver(const SolverContext<dim, degree> &_solver_context,
                                     Types::Index _solve_priority = 0)
    : SequentialSolver<dim, degree, number>(_solver_context,
                                            FieldSolveType::NonexplicitAuxiliary,
                                            _solve_priority) {};

  /**
   * @brief Destructor.
   */
  ~SequentialAuxiliarySolver() override = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialAuxiliarySolver(const SequentialAuxiliarySolver &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialAuxiliarySolver &
  operator=(const SequentialAuxiliarySolver &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialAuxiliarySolver(SequentialAuxiliarySolver &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialAuxiliarySolver &
  operator=(SequentialAuxiliarySolver &&solver_base) noexcept = delete;

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

    // Init the explicit auxiliary solvers
    for (const auto &[index, variable] : this->get_subset_attributes())
      {
        this->init_explicit_solver(variable);
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
        // Compute the update
        this->get_system_matrix()[index]->compute_nonexplicit_auxiliary_update(
          this->get_dst_solution_subset().at(index),
          this->get_src_solution_subset().at(index));

        // Scale the update by the respective (Scalar/Vector) invm.
        this->get_dst_solution_subset().at(index).at(0)->scale(
          this->get_invm_handler().get_invm(index));

        // Update the solutions
        this->get_solution_handler().update(FieldSolveType::NonexplicitAuxiliary, index);

        // Apply constraints
        this->get_constraint_handler().get_constraint(index).distribute(
          *(this->get_solution_handler().get_solution_vector(index,
                                                             DependencyType::Normal)));

        // Update the ghosts
        Timer::start_section("Update ghosts");
        this->get_solution_handler().update_ghosts();
        Timer::end_section("Update ghosts");
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
