// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/types.h>

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
class LinearSolver : public GroupSolverBase<dim, degree, number>
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
  LinearSolver(SolveGroup                                _solve_group,
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
    // Zero out the ghosts
    Timer::start_section("Zero ghosts");
    solutions.zero_out_ghosts(relative_level);
    Timer::end_section("Zero ghosts");

    // Compute residual (rhs)
    mf_operators[relative_level].compute_operator(
      solutions.get_new_solution_full_vector(relative_level));

    // Set up linear solver
    using BlockVector = SolutionHandler<dim, number>::BlockVector;
    dealii::SolverCG<BlockVector> cg_solver(0 /* this->get_solver_control() */);
    try
      {
        cg_solver.solve(lhs_mf_operators[relative_level],
                        solutions.get_change_solution_full_vector(relative_level),
                        solutions.get_new_solution_full_vector(relative_level));
      }
    catch (...) // TODO: more specific catch
      {
        ConditionalOStreams::pout_base()
          << "Warning: linear solver did not converge as per set tolerances.\n";
      }

    // Update the solutions
    solutions.update(relative_level);

    // Apply constraints
    solutions.apply_constraints(relative_level);

    // Update the ghosts
    Timer::start_section("Update ghosts");
    solutions.update_ghosts(relative_level);
    Timer::end_section("Update ghosts");
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
