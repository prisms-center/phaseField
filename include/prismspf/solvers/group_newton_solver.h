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
public:
  using GroupSolverBase = GroupSolverBase<dim, degree, number>;
  using GroupSolverBase::mf_operators;
  using GroupSolverBase::solutions;
  using GroupSolverBase::solve_group;
  using GroupSolverBase::solver_context;

  /**
   * @brief Constructor.
   */
  NewtonSolver(SolveGroup                                _solve_group,
               const SolverContext<dim, degree, number> &_solver_context)
    : GroupSolverBase(_solve_group, _solver_context)
  {}

  /**
   * @brief Initialize the solver.
   */
  void
  init() override
  {
    for (unsigned int relative_level = 0; relative_level < mf_operators.size();
         ++relative_level)
      {
        mf_operators[relative_level] = MFOperator<dim, degree, number>(
          solver_context->pde_operator,
          PDEOperator<dim, degree, number>::compute_explicit_rhs,
          solver_context->field_attributes,
          solver_context->solution_indexer,
          relative_level,
          solve_group.dependencies_rhs);
      }
    reinit();
  }

  /**
   * @brief Reinitialize the solver.
   */
  void
  reinit() override
  {
    for (unsigned int relative_level = 0; relative_level < mf_operators.size();
         ++relative_level)
      {
        mf_operators[relative_level].initialize(solve_group,
                                                solutions.get_matrix_free(relative_level),
                                                solutions.get_global_to_block_index());
      }
  }

  /**
   * @brief Solve for a single update step.
   */
  void
  solve_level(unsigned int relative_level) override
  {
    number       step_length        = number(1.0); // TODO
    number       tolerance          = number(1.0); // TODO
    unsigned int max_iterations     = 1;           // TODO
    bool         newton_unconverged = true;
    unsigned int iter               = 0;

    // TODO: setup initial guess for solution vector. Maybe use old_solution if available.
    // Newton iteration loop.
    while (newton_unconverged && iter++ < max_iterations)
      {
        // Solve for Newton update.
        LinearSolver<dim, degree, number>::solve_level(relative_level);

        // TODO: Apply zero-constraints to 'change' vector

        // Perform Newton update.
        solutions.get_solution_full_vector(relative_level)
          .add(step_length, solutions.get_change_solution_full_vector(relative_level));

        // TODO: Apply constraints to solution vector

        // Check convergence.
        number l2_norm =
          solutions.get_change_solution_full_vector(relative_level).l2_norm();
        newton_unconverged = l2_norm > tolerance;
      }
    ConditionalOStreams::pout_verbose() << iter << " Newton iterations to converge.\n\n";
    if (iter >= max_iterations)
      {
        ConditionalOStreams::pout_base() << "Warning: nonlinear solver did not "
                                            "converge as per set tolerances.\n\n";
      }
  }

  /**
   * @brief Solve for a single update step.
   */
  void
  solve() override;

  /**
   * @brief Print information about the solver to summary.log.
   */
  void
  print() override;

private:
  /**
   * @brief Matrix free operators for each level
   */
  std::vector<MFOperator<dim, degree, number>> lhs_mf_operators;
};

PRISMS_PF_END_NAMESPACE
