// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/linear_solver_identity.h>
#include <prismspf/solvers/sequential_self_nonlinear_solver.h>
#include <prismspf/solvers/sequential_solver.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include <ostream>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
SequentialSelfNonlinearSolver<dim, degree, number>::SequentialSelfNonlinearSolver(
  const SolverContext<dim, degree, number> &_solver_context,
  Types::Index                              _solve_priority)
  : SequentialSolver<dim, degree, number>(_solver_context,
                                          FieldSolveType::NonexplicitSelfnonlinear,
                                          _solve_priority)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSelfNonlinearSolver<dim, degree, number>::init()
{
  // Call the base class init
  this->SequentialSolver<dim, degree, number>::init();

  // If the solver is empty we can just return early.
  if (this->solver_is_empty())
    {
      return;
    }

  // Init the linear solvers
  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      this->init_linear_solver(variable);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSelfNonlinearSolver<dim, degree, number>::reinit()
{
  // Call the base class reinit
  this->SequentialSolver<dim, degree, number>::reinit();

  // If the solver is empty we can just return early.
  if (this->solver_is_empty())
    {
      return;
    }

  // Reinit the linear solvers
  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      this->reinit_linear_solver(variable);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSelfNonlinearSolver<dim, degree, number>::solve()
{
  // Call the base class solve
  this->SequentialSolver<dim, degree, number>::solve();

  // If the solver is empty we can just return early.
  if (this->solver_is_empty())
    {
      return;
    }

  // Solve each field
  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      // Set the convergence bool, iteration counter, and step length
      bool         unconverged = true;
      unsigned int iteration   = 0;
      const auto   step_length =
        static_cast<number>(this->get_user_inputs()
                              .get_nonlinear_solve_parameters()
                              .get_nonlinear_solve_parameters(index)
                              .step_length);

      while (unconverged)
        {
          if (this->get_user_inputs().get_output_parameters().should_output(
                this->get_user_inputs().get_temporal_discretization().get_increment()))
            {
              ConditionalOStreams::pout_summary()
                << "Nonlinear solver step: " << iteration << "\n";
            }

          // Assume the solve is converged, unless proven otherwise
          unconverged = false;

          // Set the norm of the newton update
          number newton_update_norm = 0.0;

          // Perform the linear solve with the step length
          newton_update_norm = this->solve_linear_solver(variable, step_length);

          // Check the convergence of the nonlinear solve
          if (this->get_user_inputs().get_output_parameters().should_output(
                this->get_user_inputs().get_temporal_discretization().get_increment()))
            {
              ConditionalOStreams::pout_summary()
                << "  field: " << index << " Newton update norm: " << newton_update_norm
                << "\n"
                << std::flush;
            }

          if (newton_update_norm >
              static_cast<number>(this->get_user_inputs()
                                    .get_nonlinear_solve_parameters()
                                    .get_nonlinear_solve_parameters(index)
                                    .tolerance_value))
            {
              unconverged = true;
            }

          // Check if the maximum number of iterations has been reached
          if (iteration >= this->get_user_inputs()
                             .get_nonlinear_solve_parameters()
                             .get_nonlinear_solve_parameters(index)
                             .max_iterations)
            {
              unconverged = false;
              ConditionalOStreams::pout_base() << "Warning: nonlinear solver did not "
                                                  "converge as per set tolerances.\n\n"
                                               << std::flush;
            }

          // Update the iteration counter
          iteration++;
        }

      if (variable.get_pde_type() != PDEType::Auxiliary)
        {
          // The solve will have updated the "old solution" vector with the newton update
          // so it's technically the new solution. In order to update the solutions and
          // preserve the old states we copy the old solution to the new solution and
          // update.
          *(this->get_solution_handler().get_new_solution_vector(index)) =
            *(this->get_solution_handler().get_solution_vector(index,
                                                               DependencyType::Normal));
          this->get_solution_handler().update(this->get_field_solve_type(),
                                              this->get_solve_block(),
                                              index);
        }

      // Update the ghosts
      Timer::start_section("Update ghosts");
      this->get_solution_handler().update_ghosts();
      Timer::end_section("Update ghosts");
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSelfNonlinearSolver<dim, degree, number>::print()
{
  // Print the base class information
  this->SequentialSolver<dim, degree, number>::print();
}

#include "solvers/sequential_self_nonlinear_solver.inst"

PRISMS_PF_END_NAMESPACE
