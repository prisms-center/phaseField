// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/linear_solver_base.h>

#include <prismspf/config.h>

#include <memory>
#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
linearSolverBase<dim, degree>::linearSolverBase(
  const UserInputParameters<dim>                         &_user_inputs,
  const VariableAttributes                               &_variable_attributes,
  const MatrixfreeHandler<dim>                           &_matrix_free_handler,
  const ConstraintHandler<dim, degree>                   &_constraint_handler,
  SolutionHandler<dim>                                   &_solution_handler,
  std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator)
  : user_inputs(&_user_inputs)
  , variable_attributes(&_variable_attributes)
  , matrix_free_handler(&_matrix_free_handler)
  , constraint_handler(&_constraint_handler)
  , solution_handler(&_solution_handler)
  , field_index(_variable_attributes.get_field_index())
  , residual(_solution_handler.get_new_solution_vector(field_index))
  , newton_update(
      _solution_handler.get_solution_vector(field_index, DependencyType::Change))
  , pde_operator(std::move(_pde_operator))
  , solver_control(_user_inputs.get_linear_solve_parameters()
                     .get_linear_solve_parameters(field_index)
                     .max_iterations)
{
  // Creating map to match types
  subset_attributes.emplace(field_index, *variable_attributes);

  // Create the implementation of MatrixFreeOperator with the subset of variable
  // attributes
  system_matrix =
    std::make_unique<SystemMatrixType>(subset_attributes, pde_operator, field_index);
  update_system_matrix =
    std::make_unique<SystemMatrixType>(subset_attributes, pde_operator, field_index);

  // Create the residual subset of solution vectors and add the mapping to
  // MatrixFreeOperator
  residual_src.push_back(
    solution_handler->get_solution_vector(field_index, DependencyType::Normal));
  residual_global_to_local_solution.emplace(std::make_pair(field_index,
                                                           DependencyType::Normal),
                                            0);
  for (const auto &[variable_index, map] : variable_attributes->get_dependency_set_rhs())
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
  // MatrixFreeOperator. For this one we consider the singular src vector as the residual
  // vector above and the src subset all other dependencies vectors. This is
  // complicated, but has to do with the way deal.II's solvers (like CG) handle
  // iterative updates. For this reason, we have to pass the residual vector as
  // VectorType src and all other dependencies for the LHS as std::vector<VectorType*>
  // src_subset.

  for (const auto &[variable_index, map] : variable_attributes->get_dependency_set_lhs())
    {
      for (const auto &[dependency_type, field_type] : map)
        {
          const auto pair = std::make_pair(variable_index, dependency_type);

          if (dependency_type == DependencyType::Change)
            {
              Assert(field_index == variable_index,
                     dealii::ExcMessage("The change type should have the same type as "
                                        "the field we're solving."));
            }
          newton_update_src.push_back(
            solution_handler->get_solution_vector(variable_index, dependency_type));
          newton_update_global_to_local_solution.emplace(pair,
                                                         newton_update_src.size() - 1);
        }
    }
  Assert(
    newton_update_global_to_local_solution.size() == newton_update_src.size(),
    dealii::ExcMessage(
      "The newton update src and global to local mappings must have the same size."));
}

template <unsigned int dim, unsigned int degree>
inline void
linearSolverBase<dim, degree>::compute_solver_tolerance()
{
  tolerance = user_inputs->get_linear_solve_parameters()
                    .get_linear_solve_parameters(field_index)
                    .tolerance_type == SolverToleranceType::RelativeResidualChange
                ? user_inputs->get_linear_solve_parameters()
                      .get_linear_solve_parameters(field_index)
                      .tolerance *
                    residual->l2_norm()
                : user_inputs->get_linear_solve_parameters()
                    .get_linear_solve_parameters(field_index)
                    .tolerance;
}

INSTANTIATE_BI_TEMPLATE(linearSolverBase)

PRISMS_PF_END_NAMESPACE
