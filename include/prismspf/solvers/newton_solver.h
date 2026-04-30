// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/solver_cg.h>

#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/linear_solver.h>
#include <prismspf/solvers/mf_operator.h>
#include <prismspf/solvers/solver_base.h>

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
  using SolverBase<dim, degree, number>::solutions;
  using SolverBase<dim, degree, number>::solve_context;
  using SolverBase<dim, degree, number>::solve_block;
  using LinearSolver<dim, degree, number>::do_linear_solve;
  using LinearSolver<dim, degree, number>::normalization_value;
  using LinearSolver<dim, degree, number>::lhs_operators;
  using LinearSolver<dim, degree, number>::rhs_operators;
  using LinearSolver<dim, degree, number>::rhs_vector;

public:
  /**
   * @brief Constructor.
   */
  NewtonSolver(SolveBlock                               _solve_block,
               const SolveContext<dim, degree, number> &_solve_context)
    : LinearSolver<dim, degree, number>(_solve_block, _solve_context)
    , newton_params(
        solve_context->get_user_inputs().nonlinear_solve_parameters.newton_solvers.at(
          solve_block.id))
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
    const number newton_step_length = newton_params.step_length;
    const number newton_tolerance = newton_params.tolerance_value * normalization_value();
    unsigned int newton_max_iterations = newton_params.max_iterations;

    BlockVector<number>             &newton_residual = rhs_vector[relative_level];
    BlockVector<number>             &newton_update   = newton_updates[relative_level];
    MFOperator<dim, degree, number> &rhs_op          = rhs_operators[relative_level];
    MFOperator<dim, degree, number> &lhs_op          = lhs_operators[relative_level];
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
    number       l2_norm            = -1.0;
    int          total_lin_iters    = 0;
    while (newton_unconverged && iter < newton_max_iterations)
      {
        // Apply constraints to solution vector
        solutions.apply_constraints(relative_level);
        solutions.update_ghosts(relative_level);

        // Solve for Newton-residual (r)
        Timer::start_section("Zero ghosts");
        newton_residual.zero_out_ghost_values();
        Timer::end_section("Zero ghosts");
        rhs_op.compute_operator(newton_residual);
        newton_residual.update_ghost_values();

        // Check convergence.
        l2_norm            = newton_residual.l2_norm();
        newton_unconverged = l2_norm > newton_tolerance;
        if (!newton_unconverged || iter == newton_max_iterations)
          {
            continue;
          }

        // Solve for Newton update. (-dr/du|Du)
        total_lin_iters += do_linear_solve(newton_residual, lhs_op, newton_update);
        newton_update.update_ghost_values();

        // Zero out the ghosts
        solutions.get_solution_full_vector(relative_level).zero_out_ghost_values();

        // Perform Newton update.
        solutions.get_solution_full_vector(relative_level)
          .add(newton_step_length, newton_update);

        // Update the ghosts
        Timer::start_section("Update ghosts");
        solutions.update_ghosts(relative_level);
        Timer::end_section("Update ghosts");

        iter++;
        // Todo: implement some super simple backtracking. Something like if the residual
        // increases instead of decreasing, try again with a smaller step length.
      }
    if (solve_context->get_user_inputs().output_parameters.should_output(
          solve_context->get_simulation_timer().get_increment()))
      {
        ConditionalOStreams::pout_summary()
          << " Newton solve final residual : " << l2_norm / normalization_value()
          << " Newton steps: " << iter << " Total linear steps: " << total_lin_iters
          << "\n"
          << std::flush;
      }
    if (iter >= newton_max_iterations)
      {
        ConditionalOStreams::pout_base()
          << "[Increment " << solve_context->get_simulation_timer().get_increment()
          << "] "
          << "Warning: Newton solver did not converge before " << newton_max_iterations
          << " iterations.\n\n";
      }
  }

protected:
  std::vector<BlockVector<number>> newton_updates; //"change" term
  NonlinearSolverParameters        newton_params;
};

PRISMS_PF_END_NAMESPACE
