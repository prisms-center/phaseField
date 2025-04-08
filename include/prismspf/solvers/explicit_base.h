// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/numerics/vector_tools.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Base class for explicit solves.
 */
template <int dim, int degree>
class explicitBase
{
public:
  using SystemMatrixType = matrixFreeOperator<dim, degree, double>;

  /**
   * \brief Constructor.
   */
  explicitBase(const userInputParameters<dim> &_user_inputs,
               const matrixfreeHandler<dim>   &_matrix_free_handler,
               const invmHandler<dim, degree> &_invm_handler,
               const constraintHandler<dim>   &_constraint_handler,
               const dofHandler<dim>          &_dof_handler,
               const dealii::MappingQ1<dim>   &_mapping,
               solutionHandler<dim>           &_solution_handler,
               std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator);

  /**
   * \brief Destructor.
   */
  ~explicitBase() = default;

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
   * fieldSolveType. This function should only be used for concurrent fieldSolveTypes,
   * such as EXPLICIT, NONEXPLICIT_CO_NONLINEAR, and EXPLICIT_POSTPROCESS.
   */
  void
  compute_subset_attributes(const fieldSolveType &field_solve_type);

  /**
   * \brief Compute the shared dependency set and copy it to all eval_flag_set_RHS. Also
   * do something similar with dependency_set_RHS so that all the FEEvaluation objects are
   * initialized.
   */
  void
  compute_shared_dependencies();

  /**
   * \brief Set the initial condition according to subset_attributes.
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
  const matrixfreeHandler<dim> *matrix_free_handler;

  /**
   * \brief invm handler.
   */
  const invmHandler<dim, degree> *invm_handler;

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
   * \brief Solution handler.
   */
  solutionHandler<dim> *solution_handler;

  /**
   * \brief Subset of variable attributes.
   */
  std::map<unsigned int, variableAttributes> subset_attributes;

  /**
   * \brief PDE operator.
   */
  std::shared_ptr<const PDEOperator<dim, degree, double>> pde_operator;

  /**
   * \brief Matrix-free operator.
   */
  std::unique_ptr<SystemMatrixType> system_matrix;
};

PRISMS_PF_END_NAMESPACE