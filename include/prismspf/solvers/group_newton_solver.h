// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/solver_cg.h>

#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/group_linear_solver.h>
#include <prismspf/solvers/group_solver_base.h>
#include <prismspf/solvers/mf_operator.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolveContext;

/**
 * @brief This class handles the explicit solves of all explicit fields
 */
template <unsigned int dim, unsigned int degree, typename number>
class NewtonSolver : public LinearSolver<dim, degree, number>
{
protected:
  using GroupSolverBase<dim, degree, number>::solutions;
  using GroupSolverBase<dim, degree, number>::solve_context;
  using GroupSolverBase<dim, degree, number>::solve_group;
  using LinearSolver<dim, degree, number>::do_linear_solve;
  using LinearSolver<dim, degree, number>::lhs_operators;
  using LinearSolver<dim, degree, number>::rhs_operators;
  using LinearSolver<dim, degree, number>::rhs_vector;

public:
  /**
   * @brief Constructor.
   */
  NewtonSolver(SolveGroup                               _solve_group,
               const SolveContext<dim, degree, number> &_solve_context)
    : LinearSolver<dim, degree, number>(_solve_group, _solve_context)
  {}

  /**
   * @brief Initialize the solver.
   */
  void
  init(const std::list<DependencyMap> &all_dependeny_sets) override
  {
    LinearSolver<dim, degree, number>::init(all_dependeny_sets);
    unsigned int num_levels = solve_context->get_dof_manager().get_dof_handlers().size();
    newton_updates.resize(num_levels);
    for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
      {
        newton_updates[relative_level].reinit(
          solutions.get_solution_full_vector(relative_level));
      }
    newton_update_constraints.resize(num_levels);
    for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
      {
        newton_update_constraints[relative_level] =
          solve_context->get_constraint_manager().get_change_constraints(
            solve_group.field_indices);
      }
  }

  /**
   * @brief Reinitialize the solver.
   */
  void
  reinit() override
  {
    LinearSolver<dim, degree, number>::reinit();
    const unsigned int num_levels = newton_updates.size();
    for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
      {
        newton_updates[relative_level].reinit(
          solutions.get_solution_full_vector(relative_level));
      }
  }

  /**
   * @brief Solve for a single update step.
   */
  void
  solve_level(unsigned int relative_level) override
  {
    const NonlinearSolverParameters &params =
      solve_context->get_user_inputs().nonlinear_solve_parameters.newton_solvers.at(
        solve_group.id);
    number       newton_step_length    = params.step_length;
    number       newton_tolerance      = params.tolerance_value;
    unsigned int newton_max_iterations = params.max_iterations;

    BlockVector<number>             &newton_residual = rhs_vector[relative_level];
    BlockVector<number>             &newton_update   = newton_updates[relative_level];
    MFOperator<dim, degree, number> &rhs_op          = rhs_operators[relative_level];
    MFOperator<dim, degree, number> &lhs_op          = lhs_operators[relative_level];
    std::vector<const dealii::AffineConstraints<number> *> &change_constraints =
      newton_update_constraints[relative_level];
    // TODO: setup initial guess for solution vector. Maybe use old_solution if
    // available.
    if (!solutions.get_solution_level(relative_level).old_solutions.empty())
      {
        solutions.get_solution_full_vector(relative_level) =
          solutions.get_old_solution_full_vector(0, relative_level);
      }

    // Newton iteration loop.
    bool         newton_unconverged = true;
    unsigned int iter               = 0;
    while (newton_unconverged && iter++ < newton_max_iterations)
      {
        // Zero out the ghosts

        // Solve for Newton-residual (r)
        Timer::start_section("Zero ghosts");
        newton_residual.zero_out_ghost_values();
        Timer::end_section("Zero ghosts");
        rhs_op.compute_operator(newton_residual);
        newton_residual.update_ghost_values(); // needed?

        // Solve for Newton update. (-dr/du|Du)
        do_linear_solve(newton_residual, lhs_op, newton_update, relative_level);
        newton_update.update_ghost_values(); // needed?

        // Perform Newton update.
        solutions.get_solution_full_vector(relative_level)
          .add(newton_step_length, newton_update);

        // Apply constraints to solution vector
        solutions.apply_constraints(relative_level);

        // Update the ghosts
        Timer::start_section("Update ghosts");
        solutions.update_ghosts(relative_level);
        Timer::end_section("Update ghosts");

        // Check convergence.
        number l2_norm     = newton_residual.l2_norm();
        newton_unconverged = l2_norm > newton_tolerance;

        // Todo: implement some super simple backtracking. Something like if the residual
        // increases instead of decreasing, try again with a smaller step length.
      }
    ConditionalOStreams::pout_verbose() << iter << " Newton iterations to converge.\n\n";
    if (iter >= newton_max_iterations)
      {
        ConditionalOStreams::pout_base() << "Warning: nonlinear solver did not "
                                            "converge as per set tolerances.\n\n";
      }
  }

protected:
  std::vector<BlockVector<number>> newton_updates; //"change" term
  // TODO: consider giving ConstraintHandler a function to apply change-constraints rather
  // than keeping them here
  std::vector<std::vector<const dealii::AffineConstraints<number> *>>
    newton_update_constraints;
};

PRISMS_PF_END_NAMESPACE
