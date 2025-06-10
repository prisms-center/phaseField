// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/linear_solver_identity.h>
#include <prismspf/solvers/nonexplicit_base.h>
#include <prismspf/solvers/nonexplicit_self_nonlinear_solver.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
NonexplicitSelfnonlinearSolver<dim, degree>::NonexplicitSelfnonlinearSolver(
  const SolverContext<dim, degree> &_solver_context,
  const MGInfo<dim>                &_mg_info)
  : NonexplicitBase<dim, degree>(_solver_context)
  , mg_info(&_mg_info)
{}

template <unsigned int dim, unsigned int degree>
void
NonexplicitSelfnonlinearSolver<dim, degree>::init()
{
  this->compute_subset_attributes(FieldSolveType::NonexplicitSelfnonlinear);

  // If the subset attribute is empty return early
  if (this->get_subset_attributes().empty())
    {
      return;
    }

  this->set_initial_condition();

  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      if (this->get_user_inputs()
            .get_linear_solve_parameters()
            .get_linear_solve_parameters(index)
            .preconditioner == PreconditionerType::GMG)
        {
          gmg_solvers.emplace(
            index,
            std::make_unique<GMGSolver<dim, degree>>(this->get_user_inputs(),
                                                     variable,
                                                     this->get_matrix_free_handler(),
                                                     this->get_constraint_handler(),
                                                     this->get_triangulation_handler(),
                                                     this->get_dof_handler(),
                                                     this->get_mg_matrix_free_handler(),
                                                     this->get_solution_handler(),
                                                     this->get_pde_operator(),
                                                     this->get_pde_operator_float(),
                                                     *mg_info));
          gmg_solvers.at(index)->init();
        }
      else
        {
          identity_solvers.emplace(
            index,
            std::make_unique<IdentitySolver<dim, degree>>(this->get_user_inputs(),
                                                          variable,
                                                          this->get_matrix_free_handler(),
                                                          this->get_constraint_handler(),
                                                          this->get_solution_handler(),
                                                          this->get_pde_operator()));
          identity_solvers.at(index)->init();
        }
    }
}

template <unsigned int dim, unsigned int degree>
void
NonexplicitSelfnonlinearSolver<dim, degree>::solve()
{
  // If the subset attribute is empty return early
  if (this->get_subset_attributes().empty())
    {
      return;
    }

  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      // Skip if the field type is ImplicitTimeDependent and the current increment is 0.
      if (variable.get_pde_type() == PDEType::ImplicitTimeDependent &&
          this->get_user_inputs().get_temporal_discretization().get_current_increment() ==
            0)
        {
          continue;
        }

      bool         is_converged = true;
      unsigned int iteration    = 0;
      const auto  &step_length  = this->get_user_inputs()
                                  .get_nonlinear_solve_parameters()
                                  .get_nonlinear_solve_parameters(index)
                                  .step_length;

      while (is_converged)
        {
          is_converged = false;

          // Perform the linear solve with the step length
          if (this->get_user_inputs()
                .get_linear_solve_parameters()
                .get_linear_solve_parameters(index)
                .preconditioner == PreconditionerType::GMG)
            {
              gmg_solvers.at(index)->solve(step_length);
            }
          else
            {
              identity_solvers.at(index)->solve(step_length);
            }

          iteration++;

          if (iteration < this->get_user_inputs()
                            .get_nonlinear_solve_parameters()
                            .get_nonlinear_solve_parameters(index)
                            .max_iterations)
            {
              is_converged = true;
            }
        }

      // Update the solutions
      this->get_solution_handler().update(FieldSolveType::NonexplicitSelfnonlinear,
                                          index);

      // Update the ghosts
      this->get_solution_handler().update_ghosts();
    }
}

INSTANTIATE_BI_TEMPLATE(NonexplicitSelfnonlinearSolver)

PRISMS_PF_END_NAMESPACE
