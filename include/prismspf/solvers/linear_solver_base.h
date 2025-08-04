// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/solver_control.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolverContext;

struct VariableAttributes;

template <unsigned int dim, unsigned int degree, typename number>
class InitialCondition;

template <unsigned int dim, unsigned int degree, typename number>
class MatrixFreeOperator;

template <unsigned int dim>
class UserInputParameters;

template <unsigned int dim, typename number>
class MatrixFreeContainer;

template <unsigned int dim>
class TriangulationHandler;

template <unsigned int dim, unsigned int degree, typename number>
class InvmHandler;

template <unsigned int dim, unsigned int degree, typename number>
class ConstraintHandler;

template <unsigned int dim>
class DofHandler;

template <unsigned int dim>
class MGInfo;

template <unsigned int dim, typename number>
class SolutionHandler;

template <unsigned int dim, unsigned int degree, typename number>
class PDEOperator;

/**
 * @brief Base class that handles the assembly and linear solving of a field.
 */
template <unsigned int dim, unsigned int degree, typename number>
class LinearSolverBase
{
public:
  using SystemMatrixType = MatrixFreeOperator<dim, degree, number>;
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<number>;

  /**
   * @brief Constructor.
   */
  LinearSolverBase(const SolverContext<dim, degree, number> &_solver_context,
                   const VariableAttributes                 &_variable_attributes);

  /**
   * @brief Destructor.
   */
  virtual ~LinearSolverBase() = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  LinearSolverBase(const LinearSolverBase &solver) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  LinearSolverBase &
  operator=(const LinearSolverBase &solver) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  LinearSolverBase(LinearSolverBase &&solver) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  LinearSolverBase &
  operator=(LinearSolverBase &&solver) noexcept = delete;

  /**
   * @brief Initialize the system.
   */
  virtual void
  init() = 0;

  /**
   * @brief Reinitialize the system.
   */
  virtual void
  reinit() = 0;

  /**
   * @brief Solve the system Ax=b.
   */
  virtual void
  solve(const number &step_length = 1.0) = 0;

  /**
   * @brief Get the l2-norm of the newton update.
   */
  [[nodiscard]] number
  get_newton_update_l2_norm() const
  {
    AssertThrow(newton_update != nullptr, dealii::ExcNotInitialized());
    return newton_update->l2_norm();
  };

protected:
  /**
   * @brief Compute the solver tolerance based on the specified tolerance type.
   */
  void
  compute_solver_tolerance();

  /**
   * @brief Get the user-inputs.
   */
  [[nodiscard]] const UserInputParameters<dim> &
  get_user_inputs() const
  {
    return solver_context->get_user_inputs();
  }

  /**
   * @brief Get the variable attributes.
   */
  [[nodiscard]] const VariableAttributes &
  get_variable_attributes() const
  {
    return *variable_attributes;
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
   * @brief Get the matrix-free container.
   */
  [[nodiscard]] const MatrixFreeContainer<dim, number> &
  get_matrix_free_container() const
  {
    return solver_context->get_matrix_free_container();
  }

  /**
   * @brief Get the constraint handler.
   */
  [[nodiscard]] const ConstraintHandler<dim, degree, number> &
  get_constraint_handler() const
  {
    return solver_context->get_constraint_handler();
  }

  /**
   * @brief Get the solution handler.
   */
  [[nodiscard]] SolutionHandler<dim, number> &
  get_solution_handler() const
  {
    return solver_context->get_solution_handler();
  }

  /**
   * @brief Get the field index.
   */
  [[nodiscard]] unsigned int
  get_field_index() const
  {
    return field_index;
  }

  /**
   * @brief Get the mapping from global solution vectors to the local ones for the
   * residual solve.
   */
  [[nodiscard]] const std::vector<Types::Index> &
  get_residual_global_to_local_solution() const
  {
    return residual_global_to_local_solution;
  }

  /**
   * @brief Get the subset of fields that are necessary for the source of the residual
   * solve.
   */
  [[nodiscard]] const std::vector<VectorType *> &
  get_residual_src() const
  {
    return residual_src;
  }

  /**
   * @brief Get the residual vector.
   */
  [[nodiscard]] VectorType *
  get_residual() const
  {
    return residual;
  }

  /**
   * @brief Get the mapping from global solution vectors to the local ones for the
   * newton update.
   */
  [[nodiscard]] const std::vector<Types::Index> &
  get_newton_update_global_to_local_solution() const
  {
    return newton_update_global_to_local_solution;
  }

  /**
   * @brief Get the subset of fields that are necessary for the source of the newton
   * update.
   */
  [[nodiscard]] const std::vector<VectorType *> &
  get_newton_update_src() const
  {
    return newton_update_src;
  }

  /**
   * @brief Get the newton update vector.
   */
  [[nodiscard]] VectorType *
  get_newton_update() const
  {
    return newton_update;
  }

  /**
   * @brief Get the pde operator.
   */
  [[nodiscard]] const std::shared_ptr<const PDEOperator<dim, degree, number>> &
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
   * @brief Get the system matrix.
   */
  [[nodiscard]] const std::unique_ptr<SystemMatrixType> &
  get_system_matrix() const
  {
    return system_matrix;
  }

  /**
   * @brief Get the update system matrix.
   */
  [[nodiscard]] const std::unique_ptr<SystemMatrixType> &
  get_update_system_matrix() const
  {
    return update_system_matrix;
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
   * @brief Get the solver control.
   */
  [[nodiscard]] dealii::SolverControl &
  get_solver_control()
  {
    return solver_control;
  }

  /**
   * @brief Get the solver tolerance.
   */
  [[nodiscard]] number
  get_tolerance() const
  {
    return tolerance;
  }

  /**
   * @brief Get the element volume container.
   */
  [[nodiscard]] const ElementVolumeContainer<dim, degree, number> &
  get_element_volume_container() const
  {
    return solver_context->get_element_volume_container();
  }

  /**
   * @brief Get the multigrid info.
   */
  [[nodiscard]] const MGInfo<dim> &
  get_mg_info() const
  {
    return solver_context->get_mg_info();
  }

  /**
   * @brief Get the triangulation handler.
   */
  [[nodiscard]] const TriangulationHandler<dim> &
  get_triangulation_handler() const
  {
    return solver_context->get_triangulation_handler();
  }

private:
  /**
   * @brief Solver context.
   */
  const SolverContext<dim, degree, number> *solver_context;

  /**
   * @brief Variable attributes for field.
   */
  const VariableAttributes *variable_attributes;

  /**
   * @brief The field index we are solving.
   */
  unsigned int field_index;

  /**
   * @brief Mapping from global solution vectors to the local ones for the residual solve.
   */
  std::vector<Types::Index> residual_global_to_local_solution;

  /**
   * @brief Subset of fields that are necessary for the source of the residual solve.
   */
  std::vector<VectorType *> residual_src;

  /**
   * @brief Residual vector.
   */
  VectorType *residual;

  /**
   * @brief Mapping from global solution vectors to the local ones for the newton update.
   */
  std::vector<Types::Index> newton_update_global_to_local_solution;

  /**
   * @brief Subset of fields that are necessary for the source of the newton update.
   */
  std::vector<VectorType *> newton_update_src;

  /**
   * @brief Newton update vector.
   */
  VectorType *newton_update;

  /**
   * @brief Matrix-free operator for the residual side.
   */
  std::unique_ptr<SystemMatrixType> system_matrix;

  /**
   * @brief Matrix-free operator for the newton update side.
   */
  std::unique_ptr<SystemMatrixType> update_system_matrix;

  /**
   * @brief Subset attributes.
   */
  std::map<unsigned int, VariableAttributes> subset_attributes;

  /**
   * @brief Solver control.
   */
  dealii::SolverControl solver_control;

  /**
   * @brief Solver tolerance
   */
  number tolerance = 0.0;
};

PRISMS_PF_END_NAMESPACE
