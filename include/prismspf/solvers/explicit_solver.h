// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/timer.h>
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
class ExplicitSolver : public SolverBase<dim, degree, number>
{
  using SolverBase<dim, degree, number>::solutions;
  using SolverBase<dim, degree, number>::solve_context;
  using SolverBase<dim, degree, number>::solve_block;

public:
  /**
   * @brief Constructor.
   */
  ExplicitSolver(SolveBlock                               _solve_block,
                 const SolveContext<dim, degree, number> &_solve_context)
    : SolverBase<dim, degree, number>(_solve_block, _solve_context)
  {}

  void
  init(const std::list<DependencyMap> &all_dependeny_sets) override
  {
    SolverBase<dim, degree, number>::init(all_dependeny_sets);
    unsigned int num_levels = solve_context->get_dof_manager().get_dof_handlers().size();
    // Initialize rhs_operators
    rhs_operators.reserve(num_levels);
    for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
      {
        rhs_operators.emplace_back(solve_context->get_pde_operator(),
                                   &PDEOperatorBase<dim, degree, number>::compute_rhs,
                                   solve_context->get_field_attributes(),
                                   solve_context->get_solution_indexer(),
                                   relative_level,
                                   solve_block.dependencies_rhs,
                                   solve_context->get_simulation_timer());
        rhs_operators[relative_level].initialize(solutions);
        rhs_operators[relative_level].set_scaling_diagonal(
          true,
          solve_context->get_invm_manager().get_invm(
            solve_context->get_field_attributes(),
            solve_block.field_indices,
            relative_level));
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

    rhs_operators[relative_level].compute_operator(
      solutions.get_solution_full_vector(relative_level));

    // Apply constraints
    solutions.apply_constraints(relative_level);

    // Update the ghosts
    Timer::start_section("Update ghosts");
    solutions.update_ghosts(relative_level);
    Timer::end_section("Update ghosts");
  }

private:
  /**
   * @brief Matrix free operators for each level
   */
  std::vector<MFOperator<dim, degree, number>> rhs_operators;
};

PRISMS_PF_END_NAMESPACE
