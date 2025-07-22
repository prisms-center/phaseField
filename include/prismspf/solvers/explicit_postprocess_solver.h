// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/timer.h>

#include <prismspf/solvers/concurrent_solver.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class handles the explicit solves of all postprocessed fields
 */
template <unsigned int dim, unsigned int degree, typename number>
class ExplicitPostprocessSolver : public ConcurrentSolver<dim, degree, number>
{
public:
  /**
   * @brief Constructor.
   */
  explicit ExplicitPostprocessSolver(const SolverContext<dim, degree> &_solver_context,
                                     Types::Index _solve_priority = 0)
    : ConcurrentSolver<dim, degree, number>(_solver_context,
                                            FieldSolveType::ExplicitPostprocess,
                                            _solve_priority) {};

  /**
   * @brief Destructor.
   */
  ~ExplicitPostprocessSolver() override = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  ExplicitPostprocessSolver(const ExplicitPostprocessSolver &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  ExplicitPostprocessSolver &
  operator=(const ExplicitPostprocessSolver &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  ExplicitPostprocessSolver(ExplicitPostprocessSolver &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  ExplicitPostprocessSolver &
  operator=(ExplicitPostprocessSolver &&solver_base) noexcept = delete;

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

    // If the solver is empty we can just return early.
    if (this->solver_is_empty())
      {
        return;
      }

    // Compute the postprocessed fields
    this->get_system_matrix()->compute_postprocess_explicit_update(
      this->get_dst_solution_subset(),
      this->get_src_solution_subset());

    // Scale the update by the respective (Scalar/Vector) invm. Note that we do this with
    // the original solution set to avoid some messy mapping.
    for (auto [index, vector] : this->get_solution_handler().get_new_solution_vector())
      {
        if (this->get_subset_attributes().find(index) !=
            this->get_subset_attributes().end())
          {
            vector->scale(this->get_invm_handler().get_invm(index));
          }
      }

    // Update the solutions
    this->get_solution_handler().update(FieldSolveType::ExplicitPostprocess);

    // Update the ghosts
    Timer::start_section("Update ghosts");
    this->get_solution_handler().update_ghosts();
    Timer::end_section("Update ghosts");
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
