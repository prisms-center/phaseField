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

template <unsigned int dim>
class constraintHandler;

template <unsigned int dim>
class dofHandler;

template <unsigned int dim>
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
    const constraintHandler<dim>                           &_constraint_handler,
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
   * \brief Compute the shared dependency set and copy it to all eval_flag_set_RHS. Also
   * do something similar with dependency_set_RHS so that all the FEEvaluation objects are
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
   * \brief Print dependency_set_RHS to summary.log
   */
  void
  print();

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
  const constraintHandler<dim> *constraint_handler;

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
   * \brief Matrix-free for the newton update side.
   */
  std::map<unsigned int, std::unique_ptr<SystemMatrixType>> update_system_matrix;
};

PRISMS_PF_END_NAMESPACE