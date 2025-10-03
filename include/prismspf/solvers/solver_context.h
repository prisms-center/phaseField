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

#include <prismspf/utilities/element_volume.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class provides context for a solver with ptrs to all the relevant
 * dependencies.
 *
 * The context in this case refers to all the finite element machinery needed to solve the
 * fields. For example, it contains the triangulation handler and pde operator that
 * evaluates the user-specified PDEs.
 */
template <unsigned int dim, unsigned int degree, typename number>
class SolverContext
{
public:
  /**
   * @brief Constructor.
   */
  SolverContext(
    const UserInputParameters<dim>                         &_user_inputs,
    const MatrixFreeContainer<dim, number>                 &_matrix_free_container,
    const TriangulationHandler<dim>                        &_triangulation_handler,
    const InvmHandler<dim, degree, number>                 &_invm_handler,
    const ConstraintHandler<dim, degree, number>           &_constraint_handler,
    const DofHandler<dim>                                  &_dof_handler,
    const dealii::MappingQ1<dim>                           &_mapping,
    const ElementVolumeContainer<dim, degree, number>      &_element_volume_container,
    const MGInfo<dim>                                      &_mg_info,
    SolutionHandler<dim, number>                           &_solution_handler,
    std::shared_ptr<const PDEOperator<dim, degree, number>> _pde_operator,
    std::shared_ptr<const PDEOperator<dim, degree, float>>  _pde_operator_float)
    : user_inputs(&_user_inputs)
    , matrix_free_container(&_matrix_free_container)
    , triangulation_handler(&_triangulation_handler)
    , invm_handler(&_invm_handler)
    , constraint_handler(&_constraint_handler)
    , dof_handler(&_dof_handler)
    , mapping(&_mapping)
    , element_volume_container(&_element_volume_container)
    , mg_info(&_mg_info)
    , solution_handler(&_solution_handler)
    , pde_operator(_pde_operator)
    , pde_operator_float(_pde_operator_float) {};

  /**
   * @brief Destructor.
   */
  ~SolverContext() = default;

  /**
   * @brief Get the user-inputs.
   */
  [[nodiscard]] const UserInputParameters<dim> &
  get_user_inputs() const
  {
    Assert(user_inputs != nullptr, dealii::ExcNotInitialized());
    return *user_inputs;
  }

  /**
   * @brief Get the matrix-free object handler for non-multigrid data.
   */
  [[nodiscard]] const MatrixFreeContainer<dim, number> &
  get_matrix_free_container() const
  {
    Assert(matrix_free_container != nullptr, dealii::ExcNotInitialized());
    return *matrix_free_container;
  }

  /**
   * @brief Get the triangulation handler.
   */
  [[nodiscard]] const TriangulationHandler<dim> &
  get_triangulation_handler() const
  {
    Assert(triangulation_handler != nullptr, dealii::ExcNotInitialized());
    return *triangulation_handler;
  }

  /**
   * @brief Get the invm handler.
   */
  [[nodiscard]] const InvmHandler<dim, degree, number> &
  get_invm_handler() const
  {
    Assert(invm_handler != nullptr, dealii::ExcNotInitialized());
    return *invm_handler;
  }

  /**
   * @brief Get the constraint handler.
   */
  [[nodiscard]] const ConstraintHandler<dim, degree, number> &
  get_constraint_handler() const
  {
    Assert(constraint_handler != nullptr, dealii::ExcNotInitialized());
    return *constraint_handler;
  }

  /**
   * @brief Get the dof handler.
   */
  [[nodiscard]] const DofHandler<dim> &
  get_dof_handler() const
  {
    Assert(dof_handler != nullptr, dealii::ExcNotInitialized());
    return *dof_handler;
  }

  /**
   * @brief Get the mapping.
   */
  [[nodiscard]] const dealii::MappingQ1<dim> &
  get_mapping() const
  {
    Assert(mapping != nullptr, dealii::ExcNotInitialized());
    return *mapping;
  }

  /**
   * @brief Get the element volume container.
   */
  [[nodiscard]] const ElementVolumeContainer<dim, degree, number> &
  get_element_volume_container() const
  {
    Assert(element_volume_container != nullptr, dealii::ExcNotInitialized());
    return *element_volume_container;
  }

  /**
   * @brief Get the multigrid info.
   */
  [[nodiscard]] const MGInfo<dim> &
  get_mg_info() const
  {
    Assert(mg_info != nullptr, dealii::ExcNotInitialized());
    return *mg_info;
  }

  /**
   * @brief Get the solution handler.
   */
  [[nodiscard]] SolutionHandler<dim, number> &
  get_solution_handler() const
  {
    Assert(solution_handler != nullptr, dealii::ExcNotInitialized());
    return *solution_handler;
  }

  /**
   * @brief Get a shared pointer to the pde operator.
   */
  [[nodiscard]] const std::shared_ptr<const PDEOperator<dim, degree, number>> &
  get_pde_operator() const
  {
    Assert(pde_operator != nullptr, dealii::ExcNotInitialized());
    return pde_operator;
  }

  /**
   * @brief Get a shared pointer to the pde operator for float precision.
   */
  [[nodiscard]] const std::shared_ptr<const PDEOperator<dim, degree, float>> &
  get_pde_operator_float() const
  {
    Assert(pde_operator_float != nullptr, dealii::ExcNotInitialized());
    return pde_operator_float;
  }

private:
  /**
   * @brief User-inputs.
   */
  const UserInputParameters<dim> *user_inputs;

  /**
   * @brief Matrix-free object container.
   */
  const MatrixFreeContainer<dim, number> *matrix_free_container;

  /**
   * @brief Triangulation handler.
   */
  const TriangulationHandler<dim> *triangulation_handler;

  /**
   * @brief invm handler.
   */
  const InvmHandler<dim, degree, number> *invm_handler;

  /**
   * @brief Constraint handler.
   */
  const ConstraintHandler<dim, degree, number> *constraint_handler;

  /**
   * @brief DoF handler.
   */
  const DofHandler<dim> *dof_handler;

  /**
   * @brief Mappings to and from reference cell.
   */
  const dealii::MappingQ1<dim> *mapping;

  /**
   * @brief Element volume container.
   */
  const ElementVolumeContainer<dim, degree, number> *element_volume_container;

  /**
   * @brief Multigrid information
   */
  const MGInfo<dim> *mg_info;

  /**
   * @brief Solution handler.
   */
  SolutionHandler<dim, number> *solution_handler;

  /**
   * @brief PDE operator.
   */
  std::shared_ptr<const PDEOperator<dim, degree, number>> pde_operator;

  /**
   * @brief PDE operator for float precision.
   */
  std::shared_ptr<const PDEOperator<dim, degree, float>> pde_operator_float;
};

PRISMS_PF_END_NAMESPACE
