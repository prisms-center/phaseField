// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/invm_manager.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

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
class LinearSolver : public SolverBase<dim, degree, number>
{
protected:
  using SolverBase<dim, degree, number>::solutions;
  using SolverBase<dim, degree, number>::solve_context;
  using SolverBase<dim, degree, number>::solve_group;

public:
  /**
   * @brief Constructor.
   */
  LinearSolver(SolveGroup                               _solve_group,
               const SolveContext<dim, degree, number> &_solve_context)
    : SolverBase<dim, degree, number>(_solve_group, _solve_context)
  {}

  /**
   * @brief Initialize the solver.
   */
  void
  init(const std::list<DependencyMap> &all_dependeny_sets) override
  {
    SolverBase<dim, degree, number>::init(all_dependeny_sets);
    unsigned int num_levels = solve_context->get_dof_manager().get_dof_handlers().size();
    rhs_vector.resize(num_levels);
    for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
      {
        rhs_vector[relative_level].reinit(
          solutions.get_solution_full_vector(relative_level));
      }
    // Initialize rhs_operators
    rhs_operators.reserve(num_levels);
    for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
      {
        rhs_operators.emplace_back(solve_context->get_pde_operator(),
                                   &PDEOperatorBase<dim, degree, number>::compute_rhs,
                                   solve_context->get_field_attributes(),
                                   solve_context->get_solution_indexer(),
                                   relative_level,
                                   solve_group.dependencies_rhs,
                                   solve_context->get_simulation_timer());
        rhs_operators[relative_level].initialize(solutions);
        rhs_operators[relative_level].set_scaling_diagonal(
          true,
          solve_context->get_invm_manager().get_invm_sqrt(
            solve_context->get_field_attributes(),
            solve_group.field_indices,
            relative_level));
      }
    // Initialize lhs_operators
    lhs_operators.reserve(num_levels);
    for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
      {
        lhs_operators.emplace_back(solve_context->get_pde_operator(),
                                   &PDEOperatorBase<dim, degree, number>::compute_lhs,
                                   solve_context->get_field_attributes(),
                                   solve_context->get_solution_indexer(),
                                   relative_level,
                                   solve_group.dependencies_lhs,
                                   solve_context->get_simulation_timer());
        lhs_operators[relative_level].initialize(solutions);
        lhs_operators[relative_level].set_scaling_diagonal(
          true,
          solve_context->get_invm_manager().get_invm_sqrt(
            solve_context->get_field_attributes(),
            solve_group.field_indices,
            relative_level));
      }
  }

  /**
   * @brief Reinitialize the solver.
   */
  void
  reinit() override
  {
    SolverBase<dim, degree, number>::reinit();
    const unsigned int num_levels = rhs_vector.size();
    for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
      {
        rhs_vector[relative_level].reinit(
          solutions.get_solution_full_vector(relative_level));
      }
  }

  /**
   * @brief Solve for a single update step.
   */
  void
  solve_level(unsigned int relative_level) override
  {
    // Zero out the ghosts
    Timer::start_section("Zero ghosts");
    solutions.zero_out_ghosts(relative_level);
    Timer::end_section("Zero ghosts");

    // Set up linear solver
    rhs_operators[relative_level].compute_operator(rhs_vector[relative_level]);
    do_linear_solve(rhs_vector[relative_level],
                    lhs_operators[relative_level],
                    solutions.get_solution_full_vector(relative_level),
                    relative_level);

    // Apply constraints
    solutions.apply_constraints(relative_level);

    // Update the ghosts
    Timer::start_section("Update ghosts");
    solutions.update_ghosts(relative_level);
    Timer::end_section("Update ghosts");
  }

  void
  do_linear_solve(BlockVector<number>             &b_vector,
                  MFOperator<dim, degree, number> &lhs_operator,
                  BlockVector<number>             &x_vector,
                  unsigned int                     relative_level)
  {
    // Linear solve
    const auto &params =
      solve_context->get_user_inputs().linear_solve_parameters.linear_solvers.at(
        solve_group.id);
    try
      {
        linear_solver_control.set_max_steps(params.max_iterations);
        linear_solver_control.set_tolerance(params.tolerance);

        dealii::SolverCG<BlockVector<number>> cg_solver(linear_solver_control);
        cg_solver.solve(lhs_operator, x_vector, b_vector, dealii::PreconditionIdentity());
        if (solve_context->get_user_inputs().output_parameters.should_output(
              solve_context->get_simulation_timer().get_increment()))
          {
            ConditionalOStreams::pout_summary()
              << " Final residual (linear solve): " << linear_solver_control.last_value()
              << " Steps: " << linear_solver_control.last_step() << "\n"
              << std::flush;
          }
      }
    catch (...) // TODO: more specific catch
      {
        ConditionalOStreams::pout_base()
          << "[Increment " << solve_context->get_simulation_timer().get_increment() << "]"
          << "Warning: linear solver did not converge as per set tolerances before "
          << params.max_iterations << " iterations.\n";
      }
  }

protected:
  /**
   * @brief Matrix free operators for each level
   */
  std::vector<MFOperator<dim, degree, number>> rhs_operators;

  std::vector<MFOperator<dim, degree, number>> lhs_operators;
  std::vector<BlockVector<number>>             rhs_vector;

private:
  /**
   * @brief Solver control. Contains max iterations and tolerance.
   */
  dealii::SolverControl linear_solver_control;
};

PRISMS_PF_END_NAMESPACE
