// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/group_explicit_solver.h>
#include <prismspf/solvers/group_solver_base.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include <vector>


PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
void
ExplicitSolver<dim, degree, number>::solve_level(unsigned int relative_level)
{
  // Zero out the ghosts
  Timer::start_section("Zero ghosts");
  solutions.zero_out_ghosts(relative_level);
  Timer::end_section("Zero ghosts");

  this->get_system_matrix()->compute_explicit_update(
    solutions.get_new_solution_full_vector(relative_level),
    solutions.get_solution_full_vector(relative_level));

  solutions.solutions[0].matrix_free->cell_loop(
    &MatrixFreeOperator<dim, degree, number>::compute_local_explicit_update,
    solver_context->matrix_free_operator,
    solutions.get_new_solution_full_vector(relative_level),
    solutions.get_solution_full_vector(relative_level),
    true);

  // TODO: if it's used for all solvers, define invm in solution handler. originally
  // this was done as a loop over fields Scale the update by the respective
  // (Scalar/Vector) invm.
  solutions.get_new_solution_full_vector(relative_level).scale(block_invm);

  // Update the solutions
  solutions.update(relative_level);

  // Apply constraints
  solutions.apply_constraints(relative_level);

  // Update the ghosts
  Timer::start_section("Update ghosts");
  solutions.update_ghosts(relative_level);
  Timer::end_section("Update ghosts");
}

template <unsigned int dim, unsigned int degree, typename number>
void
ExplicitSolver<dim, degree, number>::print()
{
  // Print the base class information
  this->ConcurrentSolver<dim, degree, number>::print();
}

#include "solvers/concurrent_explicit_solver.inst"

PRISMS_PF_END_NAMESPACE