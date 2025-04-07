// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/solver_control.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * Forward declaration for user-implemented PDE class.
 */
template <int dim, int degree, typename number>
class customPDE;

/**
 * \brief Base class that handles the assembly and linear solving of a field.
 */
template <int dim, int degree>
class linearSolverBase
{
public:
  using SystemMatrixType = customPDE<dim, degree, double>;
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<double>;

  /**
   * \brief Constructor.
   */
  linearSolverBase(const userInputParameters<dim> &_user_inputs,
                   const variableAttributes       &_variable_attributes,
                   const matrixfreeHandler<dim>   &_matrix_free_handler,
                   const constraintHandler<dim>   &_constraint_handler,
                   solutionHandler<dim>           &_solution_handler);

  /**
   * \brief Destructor.
   */
  virtual ~linearSolverBase() = default;

  /**
   * \brief Initialize the system.
   */
  virtual void
  init() = 0;

  /**
   * \brief Reinitialize the system.
   */
  virtual void
  reinit() = 0;

  /**
   * \brief Solve the system Ax=b.
   */
  virtual void
  solve(const double &step_length = 1.0) = 0;

protected:
  /**
   * \brief Compute the solver tolerance based on the specified tolerance type.
   */
  void
  compute_solver_tolerance();

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> *user_inputs;

  /**
   * \brief Variable attributes for field.
   */
  const variableAttributes *variable_attributes;

  /**
   * \brief Matrix-free object handler for non-multigrid data.
   */
  const matrixfreeHandler<dim> *matrix_free_handler;

  /**
   * \brief Constraint handler.
   */
  const constraintHandler<dim> *constraint_handler;

  /**
   * \brief Solution handler.
   */
  solutionHandler<dim> *solution_handler;

  /**
   * \brief The field index we are solving.
   */
  const unsigned int field_index;

  /**
   * \brief Mapping from global solution vectors to the local ones for the residual solve.
   */
  std::unordered_map<std::pair<unsigned int, dependencyType>, unsigned int, pairHash>
    residual_global_to_local_solution;

  /**
   * \brief Subset of fields that are necessary for the source of the residual solve.
   */
  std::vector<VectorType *> residual_src;

  /**
   * \brief Residual vector.
   */
  VectorType *residual;

  /**
   * \brief Mapping from global solution vectors to the local ones for the newton update.
   */
  std::unordered_map<std::pair<unsigned int, dependencyType>, unsigned int, pairHash>
    newton_update_global_to_local_solution;

  /**
   * \brief Subset of fields that are necessary for the source of the newton update.
   */
  std::vector<VectorType *> newton_update_src;

  /**
   * \brief Newton update vector.
   */
  VectorType *newton_update;

  /**
   * \brief PDE operator for the residual side.
   */
  std::unique_ptr<SystemMatrixType> system_matrix;

  /**
   * \brief PDE operator for the newton update side.
   */
  std::unique_ptr<SystemMatrixType> update_system_matrix;

  /**
   * \brief Subset attributes.
   */
  std::map<unsigned int, variableAttributes> subset_attributes;

  /**
   * \brief Solver control.
   */
  dealii::SolverControl solver_control;

  /**
   * \brief Solver tolerance
   */
  double tolerance = 0.0;
};

template <int dim, int degree>
linearSolverBase<dim, degree>::linearSolverBase(
  const userInputParameters<dim> &_user_inputs,
  const variableAttributes       &_variable_attributes,
  const matrixfreeHandler<dim>   &_matrix_free_handler,
  const constraintHandler<dim>   &_constraint_handler,
  solutionHandler<dim>           &_solution_handler)
  : user_inputs(&_user_inputs)
  , variable_attributes(&_variable_attributes)
  , matrix_free_handler(&_matrix_free_handler)
  , constraint_handler(&_constraint_handler)
  , solution_handler(&_solution_handler)
  , field_index(_variable_attributes.field_index)
  , residual(_solution_handler.get_new_solution_vector(field_index))
  , newton_update(
      _solution_handler.get_solution_vector(field_index, dependencyType::CHANGE))
  , solver_control(
      _user_inputs.linear_solve_parameters.linear_solve.at(field_index).max_iterations)
{
  // Creating map to match types
  subset_attributes.emplace(field_index, *variable_attributes);

  // Create the implementation of customPDE with the subset of variable attributes
  system_matrix =
    std::make_unique<SystemMatrixType>(*user_inputs, field_index, subset_attributes);
  update_system_matrix =
    std::make_unique<SystemMatrixType>(*user_inputs, field_index, subset_attributes);

  // Create the residual subset of solution vectors and add the mapping to customPDE
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
  // customPDE. For this one we consider the singular src vector as the residual
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

PRISMS_PF_END_NAMESPACE