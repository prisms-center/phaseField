// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef linear_solver_base_h
#define linear_solver_base_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/solver_control.h>

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>
#include <prismspf/user_inputs/user_input_parameters.h>

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
  solve(const double step_length = 1.0) = 0;

protected:
  /**
   * \brief Compute the solver tolerance based on the specified tolerance type.
   */
  void
  compute_solver_tolerance();

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> &user_inputs;

  /**
   * \brief Variable attributes for field.
   */
  const variableAttributes &variable_attributes;

  /**
   * \brief Matrix-free object handler for non-multigrid data.
   */
  const matrixfreeHandler<dim> &matrix_free_handler;

  /**
   * \brief Constraint handler.
   */
  const constraintHandler<dim> &constraint_handler;

  /**
   * \brief Solution handler.
   */
  solutionHandler<dim> &solution_handler;

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

PRISMS_PF_END_NAMESPACE

#endif