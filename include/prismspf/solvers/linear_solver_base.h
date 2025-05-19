// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/solver_control.h>

#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

template <unsigned int dim, unsigned int degree>
class ConstraintHandler;

template <unsigned int dim, typename number>
class MatrixfreeHandler;

template <unsigned int dim>
class SolutionHandler;

template <unsigned int dim>
class TriangulationHandler;

struct VariableAttributes;

/**
 * \brief Base class that handles the assembly and linear solving of a field.
 */
template <unsigned int dim, unsigned int degree>
class linearSolverBase
{
public:
  using SystemMatrixType = MatrixFreeOperator<dim, degree, double>;
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<double>;

  /**
   * \brief Constructor.
   */
  linearSolverBase(const UserInputParameters<dim>       &_user_inputs,
                   const VariableAttributes             &_variable_attributes,
                   const MatrixfreeHandler<dim, double> &_matrix_free_handler,
                   const ConstraintHandler<dim, degree> &_constraint_handler,
                   SolutionHandler<dim>                 &_solution_handler,
                   std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator);

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
   * \brief Get the user-inputs.
   */
  [[nodiscard]] const UserInputParameters<dim> &
  get_user_inputs() const
  {
    return *user_inputs;
  }

  /**
   * \brief Get the variable attributes.
   */
  [[nodiscard]] const VariableAttributes &
  get_variable_attributes() const
  {
    return *variable_attributes;
  }

  /**
   * \brief Get the matrix-free object handler for non-multigrid data.
   */
  [[nodiscard]] const MatrixfreeHandler<dim, double> &
  get_matrix_free_handler() const
  {
    return *matrix_free_handler;
  }

  /**
   * \brief Get the constraint handler.
   */
  [[nodiscard]] const ConstraintHandler<dim, degree> &
  get_constraint_handler() const
  {
    return *constraint_handler;
  }

  /**
   * \brief Get the solution handler.
   */
  [[nodiscard]] SolutionHandler<dim> &
  get_solution_handler() const
  {
    return *solution_handler;
  }

  /**
   * \brief Get the field index.
   */
  [[nodiscard]] unsigned int
  get_field_index() const
  {
    return field_index;
  }

  /**
   * \brief Get the mapping from global solution vectors to the local ones for the
   * residual solve.
   */
  [[nodiscard]] const std::map<std::pair<unsigned int, DependencyType>, unsigned int> &
  get_residual_global_to_local_solution() const
  {
    return residual_global_to_local_solution;
  }

  /**
   * \brief Get the subset of fields that are necessary for the source of the residual
   * solve.
   */
  [[nodiscard]] const std::vector<VectorType *> &
  get_residual_src() const
  {
    return residual_src;
  }

  /**
   * \brief Get the residual vector.
   */
  [[nodiscard]] VectorType *
  get_residual() const
  {
    return residual;
  }

  /**
   * \brief Get the mapping from global solution vectors to the local ones for the
   * newton update.
   */
  [[nodiscard]] const std::map<std::pair<unsigned int, DependencyType>, unsigned int> &
  get_newton_update_global_to_local_solution() const
  {
    return newton_update_global_to_local_solution;
  }

  /**
   * \brief Get the subset of fields that are necessary for the source of the newton
   * update.
   */
  [[nodiscard]] const std::vector<VectorType *> &
  get_newton_update_src() const
  {
    return newton_update_src;
  }

  /**
   * \brief Get the newton update vector.
   */
  [[nodiscard]] VectorType *
  get_newton_update() const
  {
    return newton_update;
  }

  /**
   * \brief Get the pde operator.
   */
  [[nodiscard]] const std::shared_ptr<const PDEOperator<dim, degree, double>> &
  get_pde_operator() const
  {
    return pde_operator;
  }

  /**
   * \brief Get the system matrix.
   */
  [[nodiscard]] const std::unique_ptr<SystemMatrixType> &
  get_system_matrix() const
  {
    return system_matrix;
  }

  /**
   * \brief Get the update system matrix.
   */
  [[nodiscard]] const std::unique_ptr<SystemMatrixType> &
  get_update_system_matrix() const
  {
    return update_system_matrix;
  }

  /**
   * \brief Get the subset attributes.
   */
  [[nodiscard]] const std::map<unsigned int, VariableAttributes> &
  get_subset_attributes() const
  {
    return subset_attributes;
  }

  /**
   * \brief Get the solver control.
   */
  [[nodiscard]] dealii::SolverControl &
  get_solver_control()
  {
    return solver_control;
  }

  /**
   * \brief Get the solver tolerance.
   */
  [[nodiscard]] double
  get_tolerance() const
  {
    return tolerance;
  }

private:
  /**
   * \brief User-inputs.
   */
  const UserInputParameters<dim> *user_inputs;

  /**
   * \brief Variable attributes for field.
   */
  const VariableAttributes *variable_attributes;

  /**
   * \brief Matrix-free object handler for non-multigrid data.
   */
  const MatrixfreeHandler<dim, double> *matrix_free_handler;

  /**
   * \brief Constraint handler.
   */
  const ConstraintHandler<dim, degree> *constraint_handler;

  /**
   * \brief Solution handler.
   */
  SolutionHandler<dim> *solution_handler;

  /**
   * \brief The field index we are solving.
   */
  unsigned int field_index;

  /**
   * \brief Mapping from global solution vectors to the local ones for the residual solve.
   */
  std::map<std::pair<unsigned int, DependencyType>, unsigned int>
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
  std::map<std::pair<unsigned int, DependencyType>, unsigned int>
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
   * \brief PDE operator.
   */
  std::shared_ptr<const PDEOperator<dim, degree, double>> pde_operator;

  /**
   * \brief Matrix-free operator for the residual side.
   */
  std::unique_ptr<SystemMatrixType> system_matrix;

  /**
   * \brief Matrix-free operator for the newton update side.
   */
  std::unique_ptr<SystemMatrixType> update_system_matrix;

  /**
   * \brief Subset attributes.
   */
  std::map<unsigned int, VariableAttributes> subset_attributes;

  /**
   * \brief Solver control.
   */
  dealii::SolverControl solver_control;

  /**
   * \brief Solver tolerance
   */
  double tolerance = 0.0;
};

PRISMS_PF_END_NAMESPACE
