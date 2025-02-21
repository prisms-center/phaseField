// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/config.h>
#include <prismspf/solvers/linear_solver_base.h>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim, int degree>
linearSolverBase<dim, degree>::linearSolverBase(
  const userInputParameters<dim> &_user_inputs,
  const variableAttributes       &_variable_attributes,
  const matrixfreeHandler<dim>   &_matrix_free_handler,
  const constraintHandler<dim>   &_constraint_handler,
  solutionHandler<dim>           &_solution_handler)
  : user_inputs(_user_inputs)
  , variable_attributes(_variable_attributes)
  , matrix_free_handler(_matrix_free_handler)
  , constraint_handler(_constraint_handler)
  , solution_handler(_solution_handler)
  , field_index(_variable_attributes.field_index)
  , residual(solution_handler.new_solution_set.at(field_index))
  , newton_update(solution_handler.solution_set.at(
      std::make_pair(field_index, dependencyType::CHANGE)))
  , solver_control(
      _user_inputs.linear_solve_parameters.linear_solve.at(field_index).max_iterations)
{
  // Creating map to match types
  subset_attributes.emplace(field_index, variable_attributes);

  // Create the implementation of customPDE with the subset of variable attributes
  system_matrix =
    std::make_unique<SystemMatrixType>(user_inputs, subset_attributes, field_index);
  update_system_matrix =
    std::make_unique<SystemMatrixType>(user_inputs, subset_attributes, field_index);

  // Create the residual subset of solution vectors and add the mapping to customPDE
  residual_src.push_back(solution_handler.solution_set.at(
    std::make_pair(field_index, dependencyType::NORMAL)));
  residual_global_to_local_solution.emplace(std::make_pair(field_index,
                                                           dependencyType::NORMAL),
                                            0);
  for (const auto &[variable_index, map] : variable_attributes.dependency_set_RHS)
    {
      for (const auto &[dependency_type, field_type] : map)
        {
          const auto pair = std::make_pair(variable_index, dependency_type);

          Assert(solution_handler.solution_set.find(pair) !=
                   solution_handler.solution_set.end(),
                 dealii::ExcMessage("There is no solution vector for the given index = " +
                                    std::to_string(variable_index) +
                                    " and type = " + to_string(dependency_type)));

          Assert(solution_handler.new_solution_set.find(variable_index) !=
                   solution_handler.new_solution_set.end(),
                 dealii::ExcMessage(
                   "There is no new solution vector for the given index = " +
                   std::to_string(variable_index)));

          residual_src.push_back(solution_handler.solution_set.at(pair));
          residual_global_to_local_solution.emplace(pair, residual_src.size() - 1);
        }
    }

  // Create the newton update subset of solution vectors and add the mapping to
  // customPDE. For this one we consider the singular src vector as the residual
  // vector above and the src subset all other dependencies vectors. This is
  // complicated, but has to do with the way deal.II's solvers (like CG) handle
  // iterative updates. For this reason, we have to pass the residual vector as
  // VectorType src and all other dependencies for the LHS as std::vector<VectorType*>
  // src_subset.
  newton_update_src.push_back(solution_handler.new_solution_set.at(field_index));
  newton_update_global_to_local_solution.emplace(std::make_pair(field_index,
                                                                dependencyType::CHANGE),
                                                 0);
  for (const auto &[variable_index, map] : variable_attributes.dependency_set_LHS)
    {
      for (const auto &[dependency_type, field_type] : map)
        {
          if (dependency_type == dependencyType::CHANGE)
            {
              continue;
            }
          const auto pair = std::make_pair(variable_index, dependency_type);

          Assert(solution_handler.solution_set.find(pair) !=
                   solution_handler.solution_set.end(),
                 dealii::ExcMessage("There is no solution vector for the given index = " +
                                    std::to_string(variable_index) +
                                    " and type = " + to_string(dependency_type)));

          Assert(solution_handler.new_solution_set.find(variable_index) !=
                   solution_handler.new_solution_set.end(),
                 dealii::ExcMessage(
                   "There is no new solution vector for the given index = " +
                   std::to_string(variable_index)));

          newton_update_src.push_back(solution_handler.solution_set.at(pair));
          newton_update_global_to_local_solution.emplace(pair,
                                                         newton_update_src.size() - 1);
        }
    }
}

template <int dim, int degree>
inline void
linearSolverBase<dim, degree>::compute_solver_tolerance()
{
  tolerance =
    user_inputs.linear_solve_parameters.linear_solve.at(field_index).tolerance_type ==
        solverToleranceType::RELATIVE_RESIDUAL_CHANGE
      ? user_inputs.linear_solve_parameters.linear_solve.at(field_index).tolerance *
          residual->l2_norm()
      : user_inputs.linear_solve_parameters.linear_solve.at(field_index).tolerance;
}

INSTANTIATE_BI_TEMPLATE(linearSolverBase)

PRISMS_PF_END_NAMESPACE