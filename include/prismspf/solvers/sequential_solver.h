// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/solver_base.h>

#include <prismspf/config.h>

#include <map>
#include <memory>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolverContext;

template <unsigned int dim, unsigned int degree, typename number>
class GMGSolver;

template <unsigned int dim, unsigned int degree, typename number>
class IdentitySolver;

struct VariableAttributes;

template <unsigned int dim, unsigned int degree, typename number>
class SequentialSolver : public SolverBase<dim, degree, number>
{
public:
  /**
   * @brief Constructor.
   */
  SequentialSolver(const SolverContext<dim, degree, number> &_solver_context,
                   const FieldSolveType                     &_field_solve_type,
                   Types::Index                              _solve_priority = 0);

  /**
   * @brief Destructor.
   */
  ~SequentialSolver() override = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialSolver(const SequentialSolver &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialSolver &
  operator=(const SequentialSolver &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialSolver(SequentialSolver &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialSolver &
  operator=(SequentialSolver &&solver_base) noexcept = delete;

  /**
   * @brief Initialize the solver.
   */
  void
  init() override;

  /**
   * @brief Reinitialize the solver.
   */
  void
  reinit() override;

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

  /**
   * @brief Init a linear solver object of a given VariableAttributes.
   *
   * @param[in] variable The VariableAttributes
   */
  void
  init_linear_solver(const VariableAttributes &variable);

  /**
   * @brief Init a explicit solver objects of a given VariableAttributes.
   *
   * @param[in] variable The VariableAttributes
   */
  void
  init_explicit_solver(const VariableAttributes &variable);

  /**
   * @brief Reinit a linear solver object of a given VariableAttributes.
   *
   * @param[in] variable The VariableAttributes
   */
  void
  reinit_linear_solver(const VariableAttributes &variable);

  /**
   * @brief Reinit a explicit solver objects of a given VariableAttributes.
   *
   * @param[in] variable The VariableAttributes
   */
  void
  reinit_explicit_solver(const VariableAttributes &variable);

  /**
   * @brief Solve the explicit solver objects of a given VariableAttributes.
   *
   * @param[in] variable The VariableAttributes
   */
  void
  solve_explicit_solver(const VariableAttributes &variable);

  /**
   * @brief Solve the linear solver objects of a given VariableAttributes.
   *
   * @param[in] variable The VariableAttributes
   */
  void
  solve_linear_solver(const VariableAttributes &variable);

  /**
   * @brief Solve the linear solver objects of a given VariableAttributes.
   *
   * This function is overload specialized for co-nonlinear solves to use a given step
   * length and return a norm of the newton update. Additionally, it doesn't update the
   * solution.
   *
   * @param[in] variable The VariableAttributes
   * @param[in] step_length The step length of the linear solve. This is only used for
   * nonlinear solves when we don't want to use the entire solution.
   */
  number
  solve_linear_solver(const VariableAttributes &variable, const number &step_length);

  /**
   * @brief Get the matrix-free operator for the residual side.
   */
  [[nodiscard]] std::map<
    Types::Index,
    std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>> &
  get_system_matrix()
  {
    return system_matrix;
  }

  /**
   * @brief Get the matrix-free operator for the newton update side.
   */
  [[nodiscard]] std::map<
    Types::Index,
    std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>> &
  get_update_system_matrix()
  {
    return update_system_matrix;
  }

private:
  /**
   * @brief Matrix-free operator for the residual side.
   */
  std::map<Types::Index,
           std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>>
    system_matrix;

  /**
   * @brief Matrix-free operator for the newton update side.
   */
  std::map<Types::Index,
           std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>>
    update_system_matrix;

  /**
   * @brief Mapping from global solution vectors to the local ones
   */
  std::map<Types::Index, std::vector<Types::Index>> global_to_local_solution;

  /**
   * @brief Subset of solutions fields that are necessary for concurrent solves.
   */
  std::map<Types::Index,
           std::vector<typename SolverBase<dim, degree, number>::VectorType *>>
    solution_subset;

  /**
   * @brief Subset of new solutions fields that are necessary for concurrent solves.
   */
  std::map<Types::Index,
           std::vector<typename SolverBase<dim, degree, number>::VectorType *>>
    new_solution_subset;

  /**
   * @brief List of subset attributes.
   */
  std::vector<std::map<Types::Index, VariableAttributes>> subset_attributes_list;

  /**
   * @brief Map of identity linear solvers
   */
  std::map<Types::Index, std::unique_ptr<IdentitySolver<dim, degree, number>>>
    identity_solvers;

  /**
   * @brief Map of geometric multigrid linear solvers
   */
  std::map<Types::Index, std::unique_ptr<GMGSolver<dim, degree, number>>> gmg_solvers;
};

PRISMS_PF_END_NAMESPACE
