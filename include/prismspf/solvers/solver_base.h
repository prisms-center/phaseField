// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/numerics/vector_tools.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <map>
#include <memory>

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

template <unsigned int dim, unsigned int degree, typename number>
class ElementVolumeContainer;

template <unsigned int dim, typename number>
class SolutionHandler;

template <unsigned int dim, unsigned int degree, typename number>
class PDEOperator;

template <unsigned int dim, unsigned int degree, typename number>
class SolverBase
{
public:
  using SystemMatrixType = MatrixFreeOperator<dim, degree, number>;
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<number>;

  /**
   * @brief Constructor.
   */
  SolverBase(const SolverContext<dim, degree, number> &_solver_context,
             const FieldSolveType                     &_field_solve_type,
             Types::Index                              _solve_priority = 0);

  /**
   * @brief Destructor.
   */
  virtual ~SolverBase() = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  SolverBase(const SolverBase &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  SolverBase &
  operator=(const SolverBase &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  SolverBase(SolverBase &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  SolverBase &
  operator=(SolverBase &&solver_base) noexcept = delete;

  /**
   * @brief Initialize the solver.
   */
  virtual void
  init();

  /**
   * @brief Reinitialize the solver.
   */
  virtual void
  reinit();

  /**
   * @brief Solve for a single update step.
   */
  virtual void
  solve();

  /**
   * @brief Print information about the solver to summary.log.
   */
  virtual void
  print();

  /**
   * @brief Whether the subset attributes is empty.
   *
   * This function is used so we can return early.
   */
  [[nodiscard]] bool
  solver_is_empty() const;

  /**
   * @brief Compute the subset of VariableAttributes that belongs to a given
   * FieldSolveType and solver order.
   *
   * This function creates and returns a map of the VariablesAttributes that belong to a
   * FieldSolveType and solve order.
   */
  [[nodiscard]] std::map<Types::Index, VariableAttributes>
  compute_subset_attributes(const FieldSolveType &_field_solve_type,
                            Types::Index          _solve_priority) const;

  /**
   * @brief Compute and update the subset of VariableAttributes that belongs to a given
   * FieldSolveType and solver order.
   *
   * This function creates a map of the VariablesAttributes that belong to a
   * FieldSolveType and solve order. The map can be accessed with get_subset_attributes.
   */
  void
  update_subset_attributes(const FieldSolveType &_field_solve_type,
                           Types::Index          _solve_priority);

  /**
   * @brief Set the initial condition according to subset_attributes.
   *
   * This only sets the initial conditions for ExplicitTimeDependent,
   * ImplicitTimeDependent, TimeIndependent, and Constant fields.
   */
  void
  set_initial_condition();

  /**
   * @brief Get the user-inputs.
   */
  [[nodiscard]] const UserInputParameters<dim> &
  get_user_inputs() const
  {
    return solver_context->get_user_inputs();
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
  [[nodiscard]] const InvmHandler<dim, degree, number> &
  get_invm_handler() const
  {
    return solver_context->get_invm_handler();
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
   * @brief Get the multigrid info.
   */
  [[nodiscard]] const MGInfo<dim> &
  get_mg_info() const
  {
    return solver_context->get_mg_info();
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
   * @brief Get the element volume container.
   */
  [[nodiscard]] const ElementVolumeContainer<dim, degree, number> &
  get_element_volume_container() const
  {
    return solver_context->get_element_volume_container();
  }

  /**
   * @brief Get the subset attributes.
   */
  [[nodiscard]] const std::map<Types::Index, VariableAttributes> &
  get_subset_attributes() const
  {
    return subset_attributes;
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
   * @brief Get the field solve type.
   */
  [[nodiscard]] const FieldSolveType &
  get_field_solve_type() const
  {
    Assert(field_solve_type != Numbers::invalid_field_solve_type,
           dealii::ExcMessage("The field solve type is invalid."));
    return field_solve_type;
  }

  /**
   * @brief Get the solve block.
   */
  [[nodiscard]] Types::Index
  get_solve_block() const
  {
    return solve_priority;
  }

  /**
   * @brief Get the solver context.
   */
  [[nodiscard]] const SolverContext<dim, degree, number> &
  get_solver_context() const
  {
    Assert(solver_context != nullptr, dealii::ExcNotInitialized());
    return *solver_context;
  }

private:
  /**
   * @brief Solver context.
   */
  const std::shared_ptr<SolverContext<dim, degree, number>> solver_context;

  /**
   * @brief Field solve type.
   */
  const FieldSolveType field_solve_type = Numbers::invalid_field_solve_type;

  /**
   * @brief Solve priority.
   */
  const Types::Index solve_priority = Numbers::invalid_index;

  /**
   * @brief Subset of variable attributes.
   */
  std::map<Types::Index, VariableAttributes> subset_attributes;
};

PRISMS_PF_END_NAMESPACE
