// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/numerics/vector_tools.h>

#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class userInputParameters;

template <unsigned int dim, unsigned int degree>
class constraintHandler;

template <unsigned int dim>
class dofHandler;

template <unsigned int dim, unsigned int degree>
class initialCondition;

template <unsigned int dim, unsigned int degree, typename number>
class invmHandler;

template <unsigned int dim, typename number>
class matrixfreeHandler;

template <unsigned int dim>
class solutionHandler;

template <unsigned int dim>
class triangulationHandler;

struct variableAttributes;

/**
 * \brief Base class for nonexplicit solves.
 */
template <unsigned int dim, unsigned int degree>
class nonexplicitBase
{
public:
  using SystemMatrixType = matrixFreeOperator<dim, degree, double>;

  /**
   * \brief Constructor.
   */
  nonexplicitBase(
    const userInputParameters<dim>                         &_user_inputs,
    const matrixfreeHandler<dim, double>                   &_matrix_free_handler,
    const triangulationHandler<dim>                        &_triangulation_handler,
    const invmHandler<dim, degree, double>                 &_invm_handler,
    const constraintHandler<dim, degree>                   &_constraint_handler,
    const dofHandler<dim>                                  &_dof_handler,
    const dealii::MappingQ1<dim>                           &_mapping,
    dealii::MGLevelObject<matrixfreeHandler<dim, float>>   &_mg_matrix_free_handler,
    solutionHandler<dim>                                   &_solution_handler,
    std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator);

  /**
   * \brief Destructor.
   */
  virtual ~nonexplicitBase() = default;

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
   * \brief Compute the subset of variableAttributes that belongs to a given
   * fieldSolveType. This function should only be used for nonexplicit fieldSolveTypes,
   * such as NONEXPLICIT_LINEAR, NONEXPLICIT_SELF_NONLINEAR, NONEXPLICIT_AUXILIARY, and
   * NONEXPLICIT_CO_NONLINEAR.
   */
  void
  compute_subset_attributes(const fieldSolveType &field_solve_type);

  /**
   * \brief Compute the shared dependency set and copy it to all eval_flag_set_rhs. Also
   * do something similar with dependency_set_rhs so that all the FEEvaluation objects are
   * initialized. This should only be called for concurrent nonexplicit fieldSolveTypes
   * like NONEXPLICIT_CO_NONLINEAR.
   */
  void
  compute_shared_dependencies();

  /**
   * \brief Set the initial condition according to subset_attributes. This only applies
   * for PDEType IMPLICIT_TIME_DEPENDENT fields.
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
  [[nodiscard]] const userInputParameters<dim> &
  get_user_inputs() const
  {
    return *user_inputs;
  }

  /**
   * \brief Get the matrix-free object handler for non-multigrid data.
   */
  [[nodiscard]] const matrixfreeHandler<dim, double> &
  get_matrix_free_handler() const
  {
    return *matrix_free_handler;
  }

  /**
   * \brief Get the triangulation handler.
   */
  [[nodiscard]] const triangulationHandler<dim> &
  get_triangulation_handler() const
  {
    return *triangulation_handler;
  }

  /**
   * \brief Get the invm handler.
   */
  [[nodiscard]] const invmHandler<dim, degree, double> &
  get_invm_handler() const
  {
    return *invm_handler;
  }

  /**
   * \brief Get the constraint handler.
   */
  [[nodiscard]] const constraintHandler<dim, degree> &
  get_constraint_handler() const
  {
    return *constraint_handler;
  }

  /**
   * \brief Get the dof handler.
   */
  [[nodiscard]] const dofHandler<dim> &
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
  [[nodiscard]] dealii::MGLevelObject<matrixfreeHandler<dim, float>> &
  get_mg_matrix_free_handler()
  {
    return *mg_matrix_free_handler;
  }

  /**
   * \brief Get the solution handler.
   */
  [[nodiscard]] solutionHandler<dim> &
  get_solution_handler() const
  {
    return *solution_handler;
  }

  /**
   * \brief Get the subset attributes.
   */
  [[nodiscard]] const std::map<unsigned int, variableAttributes> &
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
  const userInputParameters<dim> *user_inputs;

  /**
   * \brief Matrix-free object handler for non-multigrid data.
   */
  const matrixfreeHandler<dim, double> *matrix_free_handler;

  /**
   * \brief Triangulation handler.
   */
  const triangulationHandler<dim> *triangulation_handler;

  /**
   * \brief invm handler.
   */
  const invmHandler<dim, degree, double> *invm_handler;

  /**
   * \brief Constraint handler.
   */
  const constraintHandler<dim, degree> *constraint_handler;

  /**
   * \brief DoF handler.
   */
  const dofHandler<dim> *dof_handler;

  /**
   * \brief Mappings to and from reference cell.
   */
  const dealii::MappingQ1<dim> *mapping;

  /**
   * \brief Matrix-free object handler for multigrid data.
   */
  dealii::MGLevelObject<matrixfreeHandler<dim, float>> *mg_matrix_free_handler;

  /**
   * \brief Solution handler.
   */
  solutionHandler<dim> *solution_handler;

  /**
   * \brief Subset of variable attributes for fields.
   */
  std::map<unsigned int, variableAttributes> subset_attributes;

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