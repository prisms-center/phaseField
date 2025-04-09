#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/linear_solver_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim, int degree>
linearSolverBase<dim, degree>::linearSolverBase(
  const userInputParameters<dim>                         &_user_inputs,
  const variableAttributes                               &_variable_attributes,
  const matrixfreeHandler<dim>                           &_matrix_free_handler,
  const constraintHandler<dim>                           &_constraint_handler,
  solutionHandler<dim>                                   &_solution_handler,
  std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator)
  : user_inputs(&_user_inputs)
  , variable_attributes(&_variable_attributes)
  , matrix_free_handler(&_matrix_free_handler)
  , constraint_handler(&_constraint_handler)
  , solution_handler(&_solution_handler)
  , field_index(_variable_attributes.field_index)
  , residual(_solution_handler.get_new_solution_vector(field_index))
  , newton_update(
      _solution_handler.get_solution_vector(field_index, dependencyType::CHANGE))
  , pde_operator(_pde_operator)
  , solver_control(
      _user_inputs.linear_solve_parameters.linear_solve.at(field_index).max_iterations)
{
  // Creating map to match types
  subset_attributes.emplace(field_index, *variable_attributes);

  // Create the implementation of matrixFreeOperator with the subset of variable
  // attributes
  system_matrix =
    std::make_unique<SystemMatrixType>(subset_attributes, pde_operator, field_index);
  update_system_matrix =
    std::make_unique<SystemMatrixType>(subset_attributes, pde_operator, field_index);

  // Create the residual subset of solution vectors and add the mapping to
  // matrixFreeOperator
  residual_src.push_back(
    solution_handler->get_solution_vector(field_index, dependencyType::NORMAL));
  residual_global_to_local_solution.emplace(std::make_pair(field_index,
                                                           dependencyType::NORMAL),
                                            0);
  for (const auto &[variable_index, map] : variable_attributes->dependency_set_RHS)
    {
      for (const auto &[dependency_type, field_type] : map)
        {
          const auto pair = std::make_pair(variable_index, dependency_type);

          residual_src.push_back(
            solution_handler->get_solution_vector(variable_index, dependency_type));
          residual_global_to_local_solution.emplace(pair, residual_src.size() - 1);
        }
    }

  // Create the newton update subset of solution vectors and add the mapping to
  // matrixFreeOperator. For this one we consider the singular src vector as the residual
  // vector above and the src subset all other dependencies vectors. This is
  // complicated, but has to do with the way deal.II's solvers (like CG) handle
  // iterative updates. For this reason, we have to pass the residual vector as
  // VectorType src and all other dependencies for the LHS as std::vector<VectorType*>
  // src_subset.
  newton_update_src.push_back(solution_handler->get_new_solution_vector(field_index));
  newton_update_global_to_local_solution.emplace(std::make_pair(field_index,
                                                                dependencyType::CHANGE),
                                                 0);
  for (const auto &[variable_index, map] : variable_attributes->dependency_set_LHS)
    {
      for (const auto &[dependency_type, field_type] : map)
        {
          if (dependency_type == dependencyType::CHANGE)
            {
              continue;
            }
          const auto pair = std::make_pair(variable_index, dependency_type);

          newton_update_src.push_back(
            solution_handler->get_solution_vector(variable_index, dependency_type));
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
    user_inputs->linear_solve_parameters.linear_solve.at(field_index).tolerance_type ==
        solverToleranceType::RELATIVE_RESIDUAL_CHANGE
      ? user_inputs->linear_solve_parameters.linear_solve.at(field_index).tolerance *
          residual->l2_norm()
      : user_inputs->linear_solve_parameters.linear_solve.at(field_index).tolerance;
}

INSTANTIATE_BI_TEMPLATE(linearSolverBase)

PRISMS_PF_END_NAMESPACE