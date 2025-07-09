// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/numerics/vector_tools.h>

#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Base class for nonexplicit solves.
 */
template <unsigned int dim, unsigned int degree>
class NonexplicitBase
{
public:
  using SystemMatrixType = MatrixFreeOperator<dim, degree, double>;

  /**
   * @brief Constructor.
   */
  explicit NonexplicitBase(const SolverContext<dim, degree> &_solver_context);

  /**
   * @brief Destructor.
   */
  virtual ~NonexplicitBase() = default;

  /**
   * @brief Initialize system.
   */
  virtual void
  init() = 0;

  /**
   * @brief Solve a single update step.
   */
  virtual void
  solve() = 0;

protected:
  /**
   * @brief Compute the subset of VariableAttributes that belongs to a given
   * FieldSolveType. This function should only be used for nonexplicit fieldSolveTypes,
   * such as NonexplicitLinear, NonexplicitSelfnonlinear, NonexplicitAuxiliary, and
   * NonexplicitCononlinear.
   */
  void
  compute_subset_attributes(const FieldSolveType &field_solve_type);

  /**
   * @brief Compute the shared dependency set and copy it to all eval_flag_set_rhs. Also
   * do something similar with dependency_set_rhs so that all the FEEvaluation objects are
   * initialized. This should only be called for concurrent nonexplicit fieldSolveTypes
   * like NonexplicitCononlinear.
   */
  void
  compute_shared_dependencies();

  /**
   * @brief Set the initial condition according to subset_attributes. This only applies
   * for PDEType ImplicitTimeDependent fields.
   */
  void
  set_initial_condition();

  /**
   * @brief Print dependency_set_rhs to summary.log
   */
  void
  print();

  /**
   * @brief Get the user-inputs.
   */
  [[nodiscard]] const UserInputParameters<dim> &
  get_user_inputs() const
  {
    return solver_context->get_user_inputs();
  }

  /**
   * @brief Get the matrix-free object handler for non-multigrid data.
   */
  [[nodiscard]] const MatrixfreeHandler<dim, double> &
  get_matrix_free_handler() const
  {
    return solver_context->get_matrix_free_handler();
  }

  /**
   * @brief Get the triangulation handler.
   */
  [[nodiscard]] const TriangulationHandler<dim> &
  get_triangulation_handler() const
  {
    return solver_context->get_triangulation_handler();
  }

  /**
   * @brief Get the invm handler.
   */
  [[nodiscard]] const InvmHandler<dim, degree, double> &
  get_invm_handler() const
  {
    return solver_context->get_invm_handler();
  }

  /**
   * @brief Get the constraint handler.
   */
  [[nodiscard]] const ConstraintHandler<dim, degree> &
  get_constraint_handler() const
  {
    return solver_context->get_constraint_handler();
  }

  /**
   * @brief Get the dof handler.
   */
  [[nodiscard]] const DofHandler<dim> &
  get_dof_handler() const
  {
    return solver_context->get_dof_handler();
  }

  /**
   * @brief Get the mapping.
   */
  [[nodiscard]] const dealii::MappingQ1<dim> &
  get_mapping() const
  {
    return solver_context->get_mapping();
  }

  /**
   * @brief Get the mg matrix-free handler.
   */
  [[nodiscard]] dealii::MGLevelObject<MatrixfreeHandler<dim, float>> &
  get_mg_matrix_free_handler()
  {
    return solver_context->get_mg_matrix_free_handler();
  }

  /**
   * @brief Get the solution handler.
   */
  [[nodiscard]] SolutionHandler<dim> &
  get_solution_handler() const
  {
    return solver_context->get_solution_handler();
  }

  /**
   * @brief Get the subset attributes.
   */
  [[nodiscard]] const std::map<unsigned int, VariableAttributes> &
  get_subset_attributes() const
  {
    return subset_attributes;
  }

  /**
   * @brief Get the pde operator.
   */
  [[nodiscard]] const std::shared_ptr<const PDEOperator<dim, degree, double>> &
  get_pde_operator() const
  {
    return solver_context->get_pde_operator();
  }

  /**
   * @brief Get the pde operator for float.
   */
  [[nodiscard]] const std::shared_ptr<const PDEOperator<dim, degree, float>> &
  get_pde_operator_float() const
  {
    return solver_context->get_pde_operator_float();
  }

  /**
   * @brief Get the matrix-free operator for the residual side.
   */
  [[nodiscard]] std::map<unsigned int, std::unique_ptr<SystemMatrixType>> &
  get_system_matrix()
  {
    return system_matrix;
  }

  /**
   * @brief Get the matrix-free operator for the newton update side.
   */
  [[nodiscard]] std::map<unsigned int, std::unique_ptr<SystemMatrixType>> &
  get_update_system_matrix()
  {
    return update_system_matrix;
  }

private:
  /**
   * @brief Solver context.
   */
  const SolverContext<dim, degree> *solver_context;

  /**
   * @brief Subset of variable attributes for fields.
   */
  std::map<unsigned int, VariableAttributes> subset_attributes;

  /**
   * @brief Matrix-free operator for the residual side.
   */
  std::map<unsigned int, std::unique_ptr<SystemMatrixType>> system_matrix;

  /**
   * @brief Matrix-free operator for the newton update side.
   */
  std::map<unsigned int, std::unique_ptr<SystemMatrixType>> update_system_matrix;
};

PRISMS_PF_END_NAMESPACE