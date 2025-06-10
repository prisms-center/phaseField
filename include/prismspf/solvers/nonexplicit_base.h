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
 * \brief Base class for nonexplicit solves.
 */
template <unsigned int dim, unsigned int degree>
class NonexplicitBase
{
public:
  using SystemMatrixType = MatrixFreeOperator<dim, degree, double>;

  /**
   * \brief Constructor.
   */
  NonexplicitBase(
    const UserInputParameters<dim>                         &_user_inputs,
    const MatrixfreeHandler<dim, double>                   &_matrix_free_handler,
    const TriangulationHandler<dim>                        &_triangulation_handler,
    const InvmHandler<dim, degree, double>                 &_invm_handler,
    const ConstraintHandler<dim, degree>                   &_constraint_handler,
    const DofHandler<dim>                                  &_dof_handler,
    const dealii::MappingQ1<dim>                           &_mapping,
    dealii::MGLevelObject<MatrixfreeHandler<dim, float>>   &_mg_matrix_free_handler,
    SolutionHandler<dim>                                   &_solution_handler,
    std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator);

  /**
   * \brief Destructor.
   */
  virtual ~NonexplicitBase() = default;

  /**
   * \brief Initialize system.
   */
  virtual void
  init() = 0;

  /**
   * \brief Solve a single update step.
   */
  virtual void
  solve() = 0;

protected:
  /**
   * \brief Compute the subset of VariableAttributes that belongs to a given
   * FieldSolveType. This function should only be used for nonexplicit fieldSolveTypes,
   * such as NonexplicitLinear, NonexplicitSelfnonlinear, NonexplicitAuxiliary, and
   * NonexplicitCononlinear.
   */
  void
  compute_subset_attributes(const FieldSolveType &field_solve_type);

  /**
   * \brief Compute the shared dependency set and copy it to all eval_flag_set_rhs. Also
   * do something similar with dependency_set_rhs so that all the FEEvaluation objects are
   * initialized. This should only be called for concurrent nonexplicit fieldSolveTypes
   * like NonexplicitCononlinear.
   */
  void
  compute_shared_dependencies();

  /**
   * \brief Set the initial condition according to subset_attributes. This only applies
   * for PDEType ImplicitTimeDependent fields.
   */
  void
  set_initial_condition();

  /**
   * \brief Print dependency_set_rhs to summary.log
   */
  void
  print();

  /**
   * \brief Get the user-inputs.
   */
  [[nodiscard]] const UserInputParameters<dim> &
  get_user_inputs() const
  {
    return *user_inputs;
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
   * \brief Get the triangulation handler.
   */
  [[nodiscard]] const TriangulationHandler<dim> &
  get_triangulation_handler() const
  {
    return *triangulation_handler;
  }

  /**
   * \brief Get the invm handler.
   */
  [[nodiscard]] const InvmHandler<dim, degree, double> &
  get_invm_handler() const
  {
    return *invm_handler;
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
   * \brief Get the dof handler.
   */
  [[nodiscard]] const DofHandler<dim> &
  get_dof_handler() const
  {
    return *dof_handler;
  }

  /**
   * \brief Get the mapping.
   */
  [[nodiscard]] const dealii::MappingQ1<dim> &
  get_mapping() const
  {
    return *mapping;
  }

  /**
   * \brief Get the mg matrix-free handler.
   */
  [[nodiscard]] dealii::MGLevelObject<MatrixfreeHandler<dim, float>> &
  get_mg_matrix_free_handler()
  {
    return *mg_matrix_free_handler;
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
   * \brief Get the subset attributes.
   */
  [[nodiscard]] const std::map<unsigned int, VariableAttributes> &
  get_subset_attributes() const
  {
    return subset_attributes;
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
   * \brief Get the matrix-free operator for the residual side.
   */
  [[nodiscard]] std::map<unsigned int, std::unique_ptr<SystemMatrixType>> &
  get_system_matrix()
  {
    return system_matrix;
  }

  /**
   * \brief Get the matrix-free operator for the newton update side.
   */
  [[nodiscard]] std::map<unsigned int, std::unique_ptr<SystemMatrixType>> &
  get_update_system_matrix()
  {
    return update_system_matrix;
  }

private:
  /**
   * \brief User-inputs.
   */
  const UserInputParameters<dim> *user_inputs;

  /**
   * \brief Matrix-free object handler for non-multigrid data.
   */
  const MatrixfreeHandler<dim, double> *matrix_free_handler;

  /**
   * \brief Triangulation handler.
   */
  const TriangulationHandler<dim> *triangulation_handler;

  /**
   * \brief invm handler.
   */
  const InvmHandler<dim, degree, double> *invm_handler;

  /**
   * \brief Constraint handler.
   */
  const ConstraintHandler<dim, degree> *constraint_handler;

  /**
   * \brief DoF handler.
   */
  const DofHandler<dim> *dof_handler;

  /**
   * \brief Mappings to and from reference cell.
   */
  const dealii::MappingQ1<dim> *mapping;

  /**
   * \brief Matrix-free object handler for multigrid data.
   */
  dealii::MGLevelObject<MatrixfreeHandler<dim, float>> *mg_matrix_free_handler;

  /**
   * \brief Solution handler.
   */
  SolutionHandler<dim> *solution_handler;

  /**
   * \brief Subset of variable attributes for fields.
   */
  std::map<unsigned int, VariableAttributes> subset_attributes;

  /**
   * \brief PDE operator.
   */
  std::shared_ptr<const PDEOperator<dim, degree, double>> pde_operator;

  /**
   * \brief Matrix-free operator for the residual side.
   */
  std::map<unsigned int, std::unique_ptr<SystemMatrixType>> system_matrix;

  /**
   * \brief Matrix-free operator for the newton update side.
   */
  std::map<unsigned int, std::unique_ptr<SystemMatrixType>> update_system_matrix;
};

PRISMS_PF_END_NAMESPACE