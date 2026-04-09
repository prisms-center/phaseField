// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/numerics/vector_tools.h>

#include <boost/geometry/core/cs.hpp>

#include <prismspf/core/dependency_extents.h>
#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/pde_operator_base.h>
#include <prismspf/core/solve_block.h>
#include <prismspf/core/system_wide.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/solve_context.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolverBase
{
public:
  /**
   * @brief Constructor.
   */
  SolverBase(SolveBlock                               _solve_block,
             const SolveContext<dim, degree, number> &_solve_context)
    : solve_block(std::move(_solve_block))
    , solve_context(&_solve_context)
    , solutions(solve_block, solve_context->get_field_attributes())
  {}

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
  init(const std::list<DependencyMap> &all_dependeny_sets)
  {
    DependencyExtents extents(solve_block.field_indices, all_dependeny_sets);
    solutions.init(solve_context->get_dof_manager(),
                   solve_context->get_constraint_manager(),
                   extents.oldest_age);

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
    solutions.reinit(solve_context->get_dof_manager(),
                     solve_context->get_constraint_manager());
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
    if (solve_context->get_simulation_timer().get_increment() == 0 &&
        solve_block.solve_timing == SolveTiming::Primary)
      {
        // Set the initial condition
        set_initial_condition();
        return;
      }
    this->solve_level(0);
    if (solve_context->get_simulation_timer().get_increment() == 0)
      {
        solutions.apply_initial_condition_for_old_fields();
      }
  }

  /**
   * @brief Update the fields.
   */
  virtual void
  update()
  {
    solutions.update();
  }

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
    for (const auto &global_index : solve_block.field_indices)
      {
        bool initialized_from_file = false;

        // First, try to find this variable in IC files
        const auto &initial_condition_parameters =
          solve_context->get_user_inputs().load_ic_parameters;

        for (const auto &initial_condition_file :
             initial_condition_parameters.get_initial_condition_files())
          {
            auto name_it =
              std::find(initial_condition_file.simulation_variable_names.begin(),
                        initial_condition_file.simulation_variable_names.end(),
                        solve_context->get_field_attributes()[global_index].name);
            if (name_it != initial_condition_file.simulation_variable_names.end())
              {
                // Found in file - read from file
                solutions.get_solution_vector(global_index).zero_out_ghost_values();
                dealii::VectorTools::interpolate(
                  SystemWide<dim, degree>::mapping,
                  solve_context->get_dof_manager().get_field_dof_handler(global_index),
                  ReadInitialCondition<dim, number>(
                    *name_it,
                    solve_context->get_field_attributes()[global_index].field_type,
                    initial_condition_file,
                    solve_context->get_user_inputs().spatial_discretization),
                  solutions.get_solution_vector(global_index));

                initialized_from_file = true;
                break; // Stop searching once found
              }
          }
        // If not found in files, use programmatic IC
        if (!initialized_from_file)
          {
            solutions.get_solution_vector(global_index).zero_out_ghost_values();
            dealii::VectorTools::interpolate(
              SystemWide<dim, degree>::mapping,
              solve_context->get_dof_manager().get_field_dof_handler(global_index),
              InitialCondition<dim, degree, number>(
                global_index,
                solve_context->get_field_attributes()[global_index].field_type,
                solve_context->get_pde_operator()),
              solutions.get_solution_vector(global_index));
          }

        solutions.apply_initial_condition_for_old_fields();
      }
  }

  /**
   * @brief Get the solution handler.
   */
  [[nodiscard]] const GroupSolutionHandler<dim, number> &
  get_solution_manager() const
  {
    return solutions;
  }

  /**
   * @brief Get the solution handler.
   */
  [[nodiscard]] GroupSolutionHandler<dim, number> &
  get_solution_manager()
  {
    return solutions;
  }

  /**
   * @brief Get the solver context.
   */
  [[nodiscard]] const SolveBlock &
  get_solve_block() const
  {
    return solve_block;
  }

protected:
  /**
   * @brief Information about the solve block this handler is responsible for.
   */
  SolveBlock solve_block;

  /**
   * @brief Solver context provides access to external information.
   */
  const SolveContext<dim, degree, number> *solve_context;

  /**
   * @brief Solution vectors for fields handled by this solver.
   */
  GroupSolutionHandler<dim, number> solutions;

  std::vector<SolverBase<dim, degree, number> *> aux_solvers;
};

PRISMS_PF_END_NAMESPACE
