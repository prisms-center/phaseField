// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/triangulation_handler.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This class provides context for a solver with ptrs to all the relevant
 * dependencies.
 */
template <unsigned int dim, unsigned int degree, typename number = double>
class SolverContext
{
public:
  /**
   * \brief Constructor.
   */
  SolverContext(
    const UserInputParameters<dim>                         &_user_inputs,
    const MatrixfreeHandler<dim, double>                   &_matrix_free_handler,
    const TriangulationHandler<dim>                        &_triangulation_handler,
    const InvmHandler<dim, degree, double>                 &_invm_handler,
    const ConstraintHandler<dim, degree>                   &_constraint_handler,
    const DofHandler<dim>                                  &_dof_handler,
    const dealii::MappingQ1<dim>                           &_mapping,
    SolutionHandler<dim>                                   &_solution_handler,
    dealii::MGLevelObject<MatrixfreeHandler<dim, float>>   &_mg_matrix_free_handler,
    std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator,
    std::shared_ptr<const PDEOperator<dim, degree, float>>  _pde_operator_float)
    : user_inputs(&_user_inputs)
    , matrix_free_handler(&_matrix_free_handler)
    , triangulation_handler(&_triangulation_handler)
    , invm_handler(&_invm_handler)
    , constraint_handler(&_constraint_handler)
    , dof_handler(&_dof_handler)
    , mapping(&_mapping)
    , solution_handler(&_solution_handler)
    , mg_matrix_free_handler(&_mg_matrix_free_handler)
    , pde_operator(std::move(_pde_operator))
    , pde_operator_float(std::move(_pde_operator_float)) {};

  /**
   * \brief Destructor.
   */
  ~SolverContext() = default;

  /**
   * \brief Get the user-inputs.
   */
  [[nodiscard]] const UserInputParameters<dim> &
  get_user_inputs() const
  {
    Assert(user_inputs != nullptr, dealii::ExcNotInitialized());
    return *user_inputs;
  }

  /**
   * \brief Get the matrix-free object handler for non-multigrid data.
   */
  [[nodiscard]] const MatrixfreeHandler<dim, double> &
  get_matrix_free_handler() const
  {
    Assert(matrix_free_handler != nullptr, dealii::ExcNotInitialized());
    return *matrix_free_handler;
  }

  /**
   * \brief Get the triangulation handler.
   */
  [[nodiscard]] const TriangulationHandler<dim> &
  get_triangulation_handler() const
  {
    Assert(triangulation_handler != nullptr, dealii::ExcNotInitialized());
    return *triangulation_handler;
  }

  /**
   * \brief Get the invm handler.
   */
  [[nodiscard]] const InvmHandler<dim, degree, double> &
  get_invm_handler() const
  {
    Assert(invm_handler != nullptr, dealii::ExcNotInitialized());
    return *invm_handler;
  }

  /**
   * \brief Get the constraint handler.
   */
  [[nodiscard]] const ConstraintHandler<dim, degree> &
  get_constraint_handler() const
  {
    Assert(constraint_handler != nullptr, dealii::ExcNotInitialized());
    return *constraint_handler;
  }

  /**
   * \brief Get the dof handler.
   */
  [[nodiscard]] const DofHandler<dim> &
  get_dof_handler() const
  {
    Assert(dof_handler != nullptr, dealii::ExcNotInitialized());
    return *dof_handler;
  }

  /**
   * \brief Get the mapping.
   */
  [[nodiscard]] const dealii::MappingQ1<dim> &
  get_mapping() const
  {
    Assert(mapping != nullptr, dealii::ExcNotInitialized());
    return *mapping;
  }

  /**
   * \brief Get the solution handler.
   */
  [[nodiscard]] SolutionHandler<dim> &
  get_solution_handler() const
  {
    Assert(solution_handler != nullptr, dealii::ExcNotInitialized());
    return *solution_handler;
  }

  /**
   * \brief Get the mg matrix-free handler.
   */
  [[nodiscard]] dealii::MGLevelObject<MatrixfreeHandler<dim, float>> &
  get_mg_matrix_free_handler() const
  {
    Assert(mg_matrix_free_handler != nullptr, dealii::ExcNotInitialized());
    return *mg_matrix_free_handler;
  }

  /**
   * \brief Get the pde operator.
   */
  [[nodiscard]] const std::shared_ptr<const PDEOperator<dim, degree, double>> &
  get_pde_operator() const
  {
    Assert(pde_operator != nullptr, dealii::ExcNotInitialized());
    return pde_operator;
  }

  /**
   * \brief Get the pde operator for float precision.
   */
  [[nodiscard]] const std::shared_ptr<const PDEOperator<dim, degree, float>> &
  get_pde_operator_float() const
  {
    Assert(pde_operator_float != nullptr, dealii::ExcNotInitialized());
    return pde_operator_float;
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
   * \brief Solution handler.
   */
  SolutionHandler<dim> *solution_handler;

  /**
   * \brief Matrix-free object handler for multigrid data.
   */
  dealii::MGLevelObject<MatrixfreeHandler<dim, float>> *mg_matrix_free_handler;

  /**
   * \brief PDE operator.
   */
  std::shared_ptr<const PDEOperator<dim, degree, double>> pde_operator;

  /**
   * \brief PDE operator for float precision.
   */
  std::shared_ptr<const PDEOperator<dim, degree, float>> pde_operator_float;
};

PRISMS_PF_END_NAMESPACE
