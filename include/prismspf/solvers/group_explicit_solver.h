// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

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
class ExplicitSolver : public GroupSolverBase<dim, degree, number>
{
  using GroupSolverBase = GroupSolverBase<dim, degree, number>;
  using GroupSolverBase::rhs_operators;
  using GroupSolverBase::solutions;
  using GroupSolverBase::solve_context;
  using GroupSolverBase::solve_group;

public:
  /**
   * @brief Constructor.
   */
  ExplicitSolver(SolveGroup                               _solve_group,
                 const SolveContext<dim, degree, number> &_solve_context)
    : GroupSolverBase(_solve_group, _solve_context)
  {}

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
      solutions.get_new_solution_full_vector(relative_level));

    // Scale by invm
    for (auto field_index : solve_group.field_indices)
      {
        solutions.get_solution_vector(field_index, relative_level)
          .scale(solve_context->get_invm_manager().get_invm(
            solve_context->get_field_attributes().at(field_index).field_type,
            relative_level));
      }

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
  // std::vector<MFOperator<dim, degree, number>> rhs_operators;
};

PRISMS_PF_END_NAMESPACE
