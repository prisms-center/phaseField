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
nonexplicitCoNonlinearSolver<dim, degree>::nonexplicitCoNonlinearSolver(
  const userInputParameters<dim>                         &_user_inputs,
  const matrixfreeHandler<dim>                           &_matrix_free_handler,
  const triangulationHandler<dim>                        &_triangulation_handler,
  const invmHandler<dim, degree>                         &_invm_handler,
  const constraintHandler<dim, degree>                   &_constraint_handler,
  const dofHandler<dim>                                  &_dof_handler,
  const dealii::MappingQ1<dim>                           &_mapping,
  dealii::MGLevelObject<matrixfreeHandler<dim, float>>   &_mg_matrix_free_handler,
  solutionHandler<dim>                                   &_solution_handler,
  std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator,
  std::shared_ptr<const PDEOperator<dim, degree, float>>  _pde_operator_float,
  const MGInfo<dim>                                      &_mg_info)
  : nonexplicitBase<dim, degree>(_user_inputs,
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
nonexplicitCoNonlinearSolver<dim, degree>::init()
{
  this->compute_subset_attributes(fieldSolveType::NONEXPLICIT_CO_NONLINEAR);

  // If the subset attribute is empty return early
  if (this->subset_attributes.empty())
    {
      return;
    }

  this->set_initial_condition();

  for (const auto &[index, variable] : this->subset_attributes)
    {
      if (variable.pde_type == PDEType::AUXILIARY)
        {
          // Creating temporary map to match types
          std::map<unsigned int, variableAttributes> temp;
          temp.emplace(index, variable);
          subset_attributes_list.push_back(temp);

          // Create the implementation of matrixFreeOperator with the subset of variable
          // attributes
          this->system_matrix[index] =
            std::make_unique<SystemMatrixType>(subset_attributes_list.back(),
                                               this->pde_operator,
                                               index);

          // Set up the user-implemented equations and create the residual vectors
          this->system_matrix.at(index)->clear();
          this->system_matrix.at(index)->initialize(
            this->matrix_free_handler->get_matrix_free());

          // Create the subset of solution vectors and add the mapping to
          // matrixFreeOperator
          new_solution_subset[index].push_back(
            this->solution_handler->get_new_solution_vector(index));
          solution_subset[index].push_back(
            this->solution_handler->get_solution_vector(index, dependencyType::NORMAL));
          global_to_local_solution[index].emplace(std::make_pair(index,
                                                                 dependencyType::NORMAL),
                                                  0);
          for (const auto &[variable_index, map] :
               subset_attributes_list.back().begin()->second.dependency_set_RHS)
            {
              for (const auto &[dependency_type, field_type] : map)
                {
                  const auto pair = std::make_pair(variable_index, dependency_type);

                  solution_subset[index].push_back(
                    this->solution_handler->get_solution_vector(variable_index,
                                                                dependency_type));
                  global_to_local_solution[index]
                    .emplace(pair, solution_subset.at(index).size() - 1);
                }
            }
          this->system_matrix.at(index)->add_global_to_local_mapping(
            global_to_local_solution.at(index));
        }
      else if (variable.pde_type == PDEType::IMPLICIT_TIME_DEPENDENT ||
               variable.pde_type == PDEType::TIME_INDEPENDENT)
        {
          if (this->user_inputs->linear_solve_parameters.linear_solve.at(index)
                .preconditioner == preconditionerType::GMG)
            {
              gmg_solvers.emplace(
                index,
                std::make_unique<GMGSolver<dim, degree>>(*this->user_inputs,
                                                         variable,
                                                         *this->matrix_free_handler,
                                                         *this->constraint_handler,
                                                         *this->triangulation_handler,
                                                         *this->dof_handler,
                                                         *this->mg_matrix_free_handler,
                                                         *this->solution_handler,
                                                         this->pde_operator,
                                                         pde_operator_float,
                                                         *mg_info));
              gmg_solvers.at(index)->init();
            }
          else
            {
              identity_solvers.emplace(
                index,
                std::make_unique<identitySolver<dim, degree>>(*this->user_inputs,
                                                              variable,
                                                              *this->matrix_free_handler,
                                                              *this->constraint_handler,
                                                              *this->solution_handler,
                                                              this->pde_operator));
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
nonexplicitCoNonlinearSolver<dim, degree>::solve()
{
  // If the subset attribute is empty return early
  if (this->subset_attributes.empty())
    {
      return;
    }

  for (const auto &[index, variable] : this->subset_attributes)
    {
      // Skip if the field type is IMPLICIT_TIME_DEPENDENT and the current increment
      // is 0.
      if (variable.pde_type == PDEType::IMPLICIT_TIME_DEPENDENT &&
          this->user_inputs->temporal_discretization.get_current_increment() == 0)
        {
          continue;
        }

      bool         is_converged = true;
      unsigned int iteration    = 0;
      const auto  &step_length =
        this->user_inputs->nonlinear_solve_parameters.nonlinear_solve.at(index)
          .step_length;

      while (is_converged)
        {
          is_converged = false;

          // Update the auxiliary fields
          if (variable.pde_type == PDEType::AUXILIARY)
            {
              // Compute the update
              this->system_matrix.at(index)->compute_nonexplicit_auxiliary_update(
                new_solution_subset.at(index),
                solution_subset.at(index));

              // Scale the update by the respective (SCALAR/VECTOR) invm.
              new_solution_subset.at(index).at(0)->scale(
                this->invm_handler->get_invm(index));

              // Update the solutions
              this->solution_handler->update(fieldSolveType::NONEXPLICIT_CO_NONLINEAR,
                                             index);

              // Apply constraints
              this->constraint_handler->get_constraint(index).distribute(
                *(this->solution_handler->get_solution_vector(index,
                                                              dependencyType::NORMAL)));
            }
          else if (variable.pde_type == PDEType::IMPLICIT_TIME_DEPENDENT ||
                   variable.pde_type == PDEType::TIME_INDEPENDENT)
            {
              // Perform the linear solve with the step length
              if (this->user_inputs->linear_solve_parameters.linear_solve.at(index)
                    .preconditioner == preconditionerType::GMG)
                {
                  gmg_solvers.at(index)->solve(step_length);
                }
              else
                {
                  identity_solvers.at(index)->solve(step_length);
                }

              iteration++;

              // TODO (landinjm): Check the convergence of the nonlinear solve somehow
              if (iteration <
                  this->user_inputs->nonlinear_solve_parameters.nonlinear_solve.at(index)
                    .max_iterations)
                {
                  is_converged = true;
                }
            }
          else
            {
              AssertThrow(false, UnreachableCode());
            }
        }
    }
}

INSTANTIATE_BI_TEMPLATE(nonexplicitCoNonlinearSolver)

PRISMS_PF_END_NAMESPACE