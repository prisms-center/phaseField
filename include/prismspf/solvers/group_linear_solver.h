// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <prismspf/core/timer.h>
#include <prismspf/core/types.h>

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
class LinearSolver : public GroupSolverBase<dim, degree, number>
{
protected:
  using GroupSolverBase<dim, degree, number>::solutions;
  using GroupSolverBase<dim, degree, number>::solve_context;
  using GroupSolverBase<dim, degree, number>::solve_group;
  using BlockVector = GroupSolutionHandler<dim, number>::BlockVector;

public:
  /**
   * @brief Constructor.
   */
  LinearSolver(SolveGroup                               _solve_group,
               const SolveContext<dim, degree, number> &_solve_context)
    : GroupSolverBase<dim, degree, number>(_solve_group, _solve_context)
  {}

  /**
   * @brief Initialize the solver.
   */
  void
  init() override
  {
    GroupSolverBase<dim, degree, number>::init();
    // Initialize rhs_operators
    for (unsigned int relative_level = 0; relative_level < rhs_operators.size();
         ++relative_level)
      {
        rhs_operators[relative_level] = MFOperator<dim, degree, number>(
          *(solve_context->get_pde_operator()),
          &PDEOperatorBase<dim, degree, number>::compute_rhs,
          solve_context->get_field_attributes(),
          solve_context->get_solution_indexer(),
          relative_level,
          solve_group.dependencies_rhs,
          solve_context->get_simulation_timer());
        rhs_operators[relative_level].initialize(solutions);
      }
    // Initialize lhs_operators
    for (unsigned int relative_level = 0; relative_level < lhs_operators.size();
         ++relative_level)
      {
        lhs_operators[relative_level] = MFOperator<dim, degree, number>(
          *(solve_context->get_pde_operator()),
          &PDEOperatorBase<dim, degree, number>::compute_lhs,
          solve_context->get_field_attributes(),
          solve_context->get_solution_indexer(),
          relative_level,
          solve_group.dependencies_lhs,
          solve_context->get_simulation_timer());
        lhs_operators[relative_level].initialize(solutions);
      }
  }

  /**
   * @brief Solve for a single update step.
   */
  void
  solve_level(unsigned int relative_level) override
  {
    // Update the solutions
    solutions.update(relative_level);

    // Zero out the ghosts
    Timer::start_section("Zero ghosts");
    solutions.zero_out_ghosts(relative_level);
    Timer::end_section("Zero ghosts");

    // Set up linear solver
    do_linear_solve(rhs_operators[relative_level],
                    solutions.get_new_solution_full_vector(relative_level),
                    lhs_operators[relative_level],
                    solutions.get_solution_full_vector(relative_level));

    // Apply constraints
    solutions.apply_constraints(relative_level);

    // Update the ghosts
    Timer::start_section("Update ghosts");
    solutions.update_ghosts(relative_level);
    Timer::end_section("Update ghosts");
  }

  void
  do_linear_solve(MFOperator<dim, degree, number> &rhs_operator,
                  BlockVector                     &b_vector,
                  MFOperator<dim, degree, number> &lhs_operator,
                  BlockVector                     &x_vector)
  {
    // Compute rhs
    rhs_operator.compute_operator(b_vector);
    // Linear solve
    try
      {
        // TODO: Make member?
        dealii::SolverCG<BlockVector> cg_solver(linear_solver_control);
        cg_solver.solve(lhs_operator, x_vector, b_vector, dealii::PreconditionIdentity());
      }
    catch (...) // TODO: more specific catch
      {
        ConditionalOStreams::pout_base()
          << "Warning: linear solver did not converge as per set tolerances.\n";
      }
  }

protected:
  /**
   * @brief Matrix free operators for each level
   */
  std::vector<MFOperator<dim, degree, number>> rhs_operators;
  std::vector<MFOperator<dim, degree, number>> lhs_operators;

private:
  /**
   * @brief Solver control. Contains max iterations and tolerance.
   */
  dealii::SolverControl linear_solver_control;
};

PRISMS_PF_END_NAMESPACE
