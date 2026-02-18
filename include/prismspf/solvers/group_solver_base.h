// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/numerics/vector_tools.h>

#include <boost/geometry/core/cs.hpp>

#include <prismspf/core/constraint_handler.h> //
#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/pde_operator_base.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/system_wide.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/solve_context.h>

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
  GroupSolverBase(SolveGroup                               _solve_group,
                  const SolveContext<dim, degree, number> &_solve_context)
    : solve_group(std::move(_solve_group))
    , solve_context(std::make_shared<SolveContext<dim, degree, number>>(_solve_context))
    , solutions(solve_group, solve_context->get_field_attributes())
  {}

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
  init()
  {
    solutions.init(SystemWide<dim, degree>::mapping,
                   solve_context->get_dof_manager(),
                   solve_context->get_constraint_manager(),
                   SystemWide<dim, degree>::quadrature);

    // Set the initial condition
    set_initial_condition();

    // Apply constraints.
    solutions.apply_constraints();
  }

  /**
   * @brief Reinitialize the solution vectors & apply constraints.
   */
  virtual void
  reinit()
  {
    // Apply constraints.
    solutions.reinit();
    solutions.apply_constraints();
  }

  /**
   * @brief Solve one level.
   */
  virtual void
  solve_level([[maybe_unused]] unsigned int relative_level = 0)
  {}

  /**
   * @brief Solve for a single update step.
   */
  virtual void
  solve()
  {
    this->solve_level(0);
  }

  /**
   * @brief Update the fields.
   */
  virtual void
  update()
  {}

  /**
   * @brief Update the ghosts.
   */
  virtual void
  update_ghosts()
  {
    solutions.update_ghosts();
  }

  /**
   * @brief Prepare for solution transfer (for AMR).
   */
  void
  prepare_for_solution_transfer()
  {
    solutions.prepare_for_solution_transfer();
  }

  /**
   * @brief Execute solution transfer (for AMR).
   */
  void
  execute_solution_transfer()
  {
    solutions.execute_solution_transfer();
  }

  /**
   * @brief Print information about the solver to summary.log.
   */
  virtual void
  print()
  {}

  /**
   * @brief Set the initial conditions.
   */
  void
  set_initial_condition()
  {
    for (const auto &global_index : solve_group.field_indices)
      {
        if (solve_context->get_user_inputs()
              .get_load_initial_condition_parameters()
              .get_read_initial_conditions_from_file())
          {
            const auto &initial_condition_parameters =
              solve_context->get_user_inputs().get_load_initial_condition_parameters();
            for (const auto &initial_condition_file :
                 initial_condition_parameters.get_initial_condition_files())
              {
                auto name_it =
                  std::find(initial_condition_file.simulation_variable_names.begin(),
                            initial_condition_file.simulation_variable_names.end(),
                            solve_context->get_field_attributes()[global_index].name);
                if (name_it != initial_condition_file.simulation_variable_names.end())
                  {
                    dealii::VectorTools::interpolate(
                      SystemWide<dim, degree>::mapping,
                      solve_context->get_dof_manager().get_dof_handler(global_index),
                      ReadInitialCondition<dim, number>(
                        *name_it,
                        solve_context->get_field_attributes()[global_index].field_type,
                        initial_condition_file,
                        solve_context->get_user_inputs().get_spatial_discretization()),
                      solutions.get_solution_vector(global_index));
                  }
              }
          }
        else
          {
            // TODO
            // dealii::VectorTools::interpolate(
            //   SystemWide<dim, degree>::mapping,
            //   solve_context->get_dof_manager().get_dof_handler(global_index),
            //   InitialCondition<dim, degree, number>(
            //     global_index,
            //     solve_context->get_field_attributes()[global_index].field_type,
            //     solve_context->get_pde_operator()),
            //   solutions.get_solution_vector(global_index));
          }
        solutions.apply_initial_condition_for_old_fields();
      }
  }

  /**
   * @brief Get the solution handler.
   */
  [[nodiscard]] GroupSolutionHandler<dim, number> &
  get_solution_manager() const;

  /**
   * @brief Get the pde operator.
   */
  [[nodiscard]] const std::shared_ptr<const PDEOperatorBase<dim, degree, number>> &
  get_pde_operator() const
  {
    return solve_context->get_pde_operator();
  }

  /**
   * @brief Get the solver context.
   */
  [[nodiscard]] const SolveContext<dim, degree, number> &
  get_solve_context() const
  {
    Assert(solve_context != nullptr, dealii::ExcNotInitialized());
    return *solve_context;
  }

protected:
  /**
   * @brief Information about the solve group this handler is responsible for.
   */
  SolveGroup solve_group;

  /**
   * @brief Solver context provides access to external information.
   */
  const std::shared_ptr<SolveContext<dim, degree, number>> solve_context;

  /**
   * @brief Solution vectors for fields handled by this solver.
   */
  GroupSolutionHandler<dim, number> solutions;

  std::vector<GroupSolverBase<dim, degree, number> *> aux_solvers;
};

PRISMS_PF_END_NAMESPACE
