// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/numerics/vector_tools.h>

#include <boost/geometry/core/cs.hpp>

#include <prismspf/core/constraint_handler.h> //
#include <prismspf/core/dof_manager.h>
#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/triangulation_manager.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/mf_operator.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class GroupSolverBase
{
public:
  /**
   * @brief Constructor.
   */
  GroupSolverBase(SolveGroup                                _solve_group,
                  const SolverContext<dim, degree, number> &_solver_context);

  /**
   * @brief Destructor.
   */
  virtual ~GroupSolverBase() = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  GroupSolverBase(const GroupSolverBase &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  GroupSolverBase &
  operator=(const GroupSolverBase &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  GroupSolverBase(GroupSolverBase &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  GroupSolverBase &
  operator=(GroupSolverBase &&solver_base) noexcept = delete;

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
   * @brief Solve one level.
   */
  virtual void
  solve_level(unsigned int relative_level = 0);

  /**
   * @brief Solve for a single update step.
   */
  virtual void
  solve();

  /**
   * @brief Update the fields. This is separate from solve because some derived solvers
   * call solve methods from other solvers.
   */
  virtual void
  update();

  /**
   * @brief Print information about the solver to summary.log.
   */
  virtual void
  print();

  /**
   * @brief Set the initial conditions.
   */
  void
  set_initial_condition();

  /**
   * @brief Get the invm handler.
   */
  [[nodiscard]] const InvmHandler<dim, degree, number> &
  get_invm_handler() const
  {
    return solver_context->get_invm_handler();
  }

  /**
   * @brief Get the solution handler.
   */
  [[nodiscard]] SolutionHandler<dim, number> &
  get_solution_manager() const;

  /**
   * @brief Get the element volume container.
   */
  [[nodiscard]] const ElementVolumeContainer<dim, degree, number> &
  get_element_volume_container() const
  {
    return solver_context->get_element_volume_container();
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
   * @brief Get the solver context.
   */
  [[nodiscard]] const SolverContext<dim, degree, number> &
  get_solver_context() const
  {
    Assert(solver_context != nullptr, dealii::ExcNotInitialized());
    return *solver_context;
  }

protected:
  /**
   * @brief List of field attributes available.
   */
  const std::vector<FieldAttributes> attributes_list;
  /**
   * @brief Information about the solve group this handler is responsible for.
   */
  SolveGroup solve_group;

  GroupSolutionHandler<dim, number>            solutions;
  std::vector<MFOperator<dim, degree, number>> mf_operators;

  std::vector<GroupSolverBase<dim, degree, number> *> aux_solvers;

  /**
   * @brief Solver context.
   */
  const std::shared_ptr<SolverContext<dim, degree, number>> solver_context;
};

PRISMS_PF_END_NAMESPACE
