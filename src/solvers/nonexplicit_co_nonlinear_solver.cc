// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/mg_level_object.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/linear_solver_identity.h>
#include <prismspf/solvers/nonexplicit_base.h>
#include <prismspf/solvers/nonexplicit_co_nonlinear_solver.h>

#include <prismspf/config.h>

#include <map>
#include <memory>
#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
NonexplicitCononlinearSolver<dim, degree>::NonexplicitCononlinearSolver(
  const UserInputParameters<dim>                         &_user_inputs,
  const MatrixfreeHandler<dim>                           &_matrix_free_handler,
  const TriangulationHandler<dim>                        &_triangulation_handler,
  const InvmHandler<dim, degree>                         &_invm_handler,
  const ConstraintHandler<dim, degree>                   &_constraint_handler,
  const DofHandler<dim>                                  &_dof_handler,
  const dealii::MappingQ1<dim>                           &_mapping,
  dealii::MGLevelObject<MatrixfreeHandler<dim, float>>   &_mg_matrix_free_handler,
  SolutionHandler<dim>                                   &_solution_handler,
  std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator,
  std::shared_ptr<const PDEOperator<dim, degree, float>>  _pde_operator_float,
  const MGInfo<dim>                                      &_mg_info)
  : NonexplicitBase<dim, degree>(_user_inputs,
                                 _matrix_free_handler,
                                 _triangulation_handler,
                                 _invm_handler,
                                 _constraint_handler,
                                 _dof_handler,
                                 _mapping,
                                 _mg_matrix_free_handler,
                                 _solution_handler,
                                 std::move(_pde_operator))
  , pde_operator_float(std::move(_pde_operator_float))
  , mg_info(&_mg_info)
{}

template <unsigned int dim, unsigned int degree>
void
NonexplicitCononlinearSolver<dim, degree>::init()
{
  this->compute_subset_attributes(FieldSolveType::NonexplicitCononlinear);

  // If the subset attribute is empty return early
  if (this->get_subset_attributes().empty())
    {
      return;
    }

  this->set_initial_condition();

  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      if (variable.get_pde_type() == PDEType::Auxiliary)
        {
          // Creating temporary map to match types
          std::map<unsigned int, VariableAttributes> temp;
          temp.emplace(index, variable);
          subset_attributes_list.push_back(temp);

          // Create the implementation of matrixFreeOperator with the subset of variable
          // attributes
          this->get_system_matrix()[index] =
            std::make_unique<SystemMatrixType>(subset_attributes_list.back(),
                                               this->get_pde_operator(),
                                               index);

          // Set up the user-implemented equations and create the residual vectors
          this->get_system_matrix().at(index)->clear();
          this->get_system_matrix().at(index)->initialize(
            this->get_matrix_free_handler().get_matrix_free());

          // Create the subset of solution vectors and add the mapping to
          // matrixFreeOperator
          new_solution_subset[index].push_back(
            this->get_solution_handler().get_new_solution_vector(index));
          solution_subset[index].push_back(
            this->get_solution_handler().get_solution_vector(index,
                                                             DependencyType::Normal));
          global_to_local_solution[index].emplace(std::make_pair(index,
                                                                 DependencyType::Normal),
                                                  0);
          for (const auto &[variable_index, map] :
               subset_attributes_list.back().begin()->second.get_dependency_set_rhs())
            {
              for (const auto &[dependency_type, field_type] : map)
                {
                  const auto pair = std::make_pair(variable_index, dependency_type);

                  solution_subset[index].push_back(
                    this->get_solution_handler().get_solution_vector(variable_index,
                                                                     dependency_type));
                  global_to_local_solution[index]
                    .emplace(pair, solution_subset.at(index).size() - 1);
                }
            }
          this->get_system_matrix().at(index)->add_global_to_local_mapping(
            global_to_local_solution.at(index));
        }
      else if (variable.get_pde_type() == PDEType::ImplicitTimeDependent ||
               variable.get_pde_type() == PDEType::TimeIndependent)
        {
          if (this->get_user_inputs()
                .get_linear_solve_parameters()
                .get_linear_solve_parameters(index)
                .preconditioner == PreconditionerType::GMG)
            {
              gmg_solvers.emplace(index,
                                  std::make_unique<GMGSolver<dim, degree>>(
                                    this->get_user_inputs(),
                                    variable,
                                    this->get_matrix_free_handler(),
                                    this->get_constraint_handler(),
                                    this->get_triangulation_handler(),
                                    this->get_dof_handler(),
                                    this->get_mg_matrix_free_handler(),
                                    this->get_solution_handler(),
                                    this->get_pde_operator(),
                                    pde_operator_float,
                                    *mg_info));
              gmg_solvers.at(index)->init();
            }
          else
            {
              identity_solvers.emplace(index,
                                       std::make_unique<IdentitySolver<dim, degree>>(
                                         this->get_user_inputs(),
                                         variable,
                                         this->get_matrix_free_handler(),
                                         this->get_constraint_handler(),
                                         this->get_solution_handler(),
                                         this->get_pde_operator()));
              identity_solvers.at(index)->init();
            }
        }
      else
        {
          AssertThrow(false, UnreachableCode());
        }
    }
}

template <unsigned int dim, unsigned int degree>
void
NonexplicitCononlinearSolver<dim, degree>::solve()
{
  // If the subset attribute is empty return early
  if (this->get_subset_attributes().empty())
    {
      return;
    }
  bool         unconverged = true;
  unsigned int iteration   = 0;

  while (unconverged)
    {
      if (this->get_user_inputs().get_output_parameters().should_output(
            this->get_user_inputs()
              .get_temporal_discretization()
              .get_current_increment()))
        {
          ConditionalOStreams::pout_summary()
            << "Nonlinear solver step: " << iteration << "\n";
        }

      // Assume the solve is converged, unless proven otherwise
      unconverged = false;

      for (const auto &[index, variable] : this->get_subset_attributes())
        {
          // Skip if the field type is IMPLICIT_TIME_DEPENDENT and the current increment
          // is 0.
          if (variable.get_pde_type() == PDEType::ImplicitTimeDependent &&
              this->get_user_inputs()
                  .get_temporal_discretization()
                  .get_current_increment() == 0)
            {
              continue;
            }

          // Get the step length
          const double step_length = this->get_user_inputs()
                                       .get_nonlinear_solve_parameters()
                                       .get_nonlinear_solve_parameters(index)
                                       .step_length;

          // Set the norm of the newton update
          double newton_update_norm = 0.0;

          // Update the auxiliary fields
          if (variable.get_pde_type() == PDEType::Auxiliary)
            {
              // Compute the update
              this->get_system_matrix().at(index)->compute_nonexplicit_auxiliary_update(
                new_solution_subset.at(index),
                solution_subset.at(index));

              // Scale the update by the respective (SCALAR/VECTOR) invm.
              new_solution_subset.at(index).at(0)->scale(
                this->get_invm_handler().get_invm(index));

              // Update the solutions
              this->get_solution_handler().update(FieldSolveType::NonexplicitCononlinear,
                                                  index);

              // Apply constraints
              this->get_constraint_handler().get_constraint(index).distribute(
                *(this->get_solution_handler()
                    .get_solution_vector(index, DependencyType::Normal)));
            }
          else if (variable.get_pde_type() == PDEType::ImplicitTimeDependent ||
                   variable.get_pde_type() == PDEType::TimeIndependent)
            {
              // Perform the linear solve with the step length
              if (this->get_user_inputs()
                    .get_linear_solve_parameters()
                    .get_linear_solve_parameters(index)
                    .preconditioner == PreconditionerType::GMG)
                {
                  gmg_solvers.at(index)->solve(step_length);
                  newton_update_norm = gmg_solvers.at(index)->get_newton_update_l2_norm();
                }
              else
                {
                  identity_solvers.at(index)->solve(step_length);
                  newton_update_norm =
                    identity_solvers.at(index)->get_newton_update_l2_norm();
                }
            }
          else
            {
              AssertThrow(false, UnreachableCode());
            }

          // Check the convergence of the nonlinear solve
          if (this->get_user_inputs().get_output_parameters().should_output(
                this->get_user_inputs()
                  .get_temporal_discretization()
                  .get_current_increment()))
            {
              ConditionalOStreams::pout_summary()
                << "  field: " << index << " Newton update norm: " << newton_update_norm
                << "\n"
                << std::flush;
            }

          if (newton_update_norm > this->get_user_inputs()
                                     .get_nonlinear_solve_parameters()
                                     .get_nonlinear_solve_parameters(index)
                                     .tolerance_value)
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
              ConditionalOStreams::pout_base()
                << "Warning: nonlinear solver did not converge as per set tolerances.\n\n"
                << std::flush;
            }
        }

      // Update the iteration counter
      iteration++;
    }

  // Update the solution
  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      this->get_solution_handler().update(FieldSolveType::NonexplicitCononlinear, index);
    }
}

INSTANTIATE_BI_TEMPLATE(NonexplicitCononlinearSolver)

PRISMS_PF_END_NAMESPACE