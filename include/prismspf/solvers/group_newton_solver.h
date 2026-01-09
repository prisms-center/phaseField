// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/lac/solver_cg.h>

#include <prismspf/core/types.h>

#include <prismspf/solvers/group_linear_solver.h>
#include <prismspf/solvers/group_solver_base.h>
#include <prismspf/solvers/mf_operator.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolverContext;

/**
 * @brief This class handles the explicit solves of all explicit fields
 */
template <unsigned int dim, unsigned int degree, typename number>
class NewtonSolver : public LinearSolver<dim, degree, number>
{
  using GroupSolverBase = GroupSolverBase<dim, degree, number>;
  using LinearSolver    = LinearSolver<dim, degree, number>;
  using GroupSolverBase::rhs_operators;
  using GroupSolverBase::solutions;
  using GroupSolverBase::solve_group;
  using GroupSolverBase::solver_context;
  using LinearSolver::do_linear_solve;
  using LinearSolver::lhs_operators;

public:
  /**
   * @brief Constructor.
   */
  NewtonSolver(SolveGroup                                _solve_group,
               const SolverContext<dim, degree, number> &_solver_context)
    : LinearSolver(_solve_group, _solver_context)
  {}

  /**
   * @brief Solve for a single update step.
   */
  void
  solve_level(unsigned int relative_level) override
  {
    number       newton_step_length    = number(1.0); // TODO
    number       newton_tolerance      = number(1.0); // TODO
    unsigned int newton_max_iterations = 1;           // TODO
    bool         newton_unconverged    = true;
    unsigned int iter                  = 0;

    // Update the solutions
    solutions.update(relative_level);

    // TODO: setup initial guess for solution vector. Maybe use old_solution if available.

    // Newton iteration loop.
    while (newton_unconverged && iter++ < newton_max_iterations)
      {
        // Zero out the ghosts
        Timer::start_section("Zero ghosts");
        solutions.zero_out_ghosts(relative_level);
        Timer::end_section("Zero ghosts");

        // Solve for Newton update.
        do_linear_solve(rhs_operators[relative_level],
                        solutions.get_new_solution_full_vector(relative_level),
                        lhs_operators[relative_level],
                        solutions.get_change_solution_full_vector(relative_level));

        // TODO: Apply zero-constraints to 'change' vector

        // Perform Newton update.
        solutions.get_solution_full_vector(relative_level)
          .add(newton_step_length,
               solutions.get_change_solution_full_vector(relative_level));

        // Apply constraints to solution vector
        solutions.apply_constraints(relative_level);

        // Update the ghosts
        Timer::start_section("Update ghosts");
        solutions.update_ghosts(relative_level);
        Timer::end_section("Update ghosts");

        // Check convergence.
        number l2_norm =
          solutions.get_change_solution_full_vector(relative_level).l2_norm();
        newton_unconverged = l2_norm > newton_tolerance;
      }
    ConditionalOStreams::pout_verbose() << iter << " Newton iterations to converge.\n\n";
    if (iter >= newton_max_iterations)
      {
        ConditionalOStreams::pout_base() << "Warning: nonlinear solver did not "
                                            "converge as per set tolerances.\n\n";
      }
  }
};

PRISMS_PF_END_NAMESPACE
