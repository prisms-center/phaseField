// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/numerics/vector_tools.h>

#include <boost/geometry/core/cs.hpp>

#include <prismspf/core/constraint_handler.h> //
#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/mf_operator.h>
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
    , solutions()
    , solve_context(std::make_shared<SolveContext<dim, degree, number>>(_solve_context))
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
    // Initialize vectors
    solutions =
      GroupSolutionHandler<dim, number>(solve_group, solve_context->field_attributes);
    solutions.init(solve_context->mapping,
                   solve_context->dof_manager,
                   solve_context->constraint_manager,
                   solve_context->quadrature);

    // Set the initial condition
    set_initial_condition();

    // Apply constraints.
    solutions.apply_constraints();

    // Initialize rhs_operators
    for (unsigned int relative_level = 0; relative_level < rhs_operators.size();
         ++relative_level)
      {
        rhs_operators[relative_level] = MFOperator<dim, degree, number>(
          solve_context->pde_operator,
          PDEOperator<dim, degree, number>::compute_explicit_rhs,
          solve_context->field_attributes,
          solve_context->solution_indexer,
          relative_level,
          solve_group.dependencies_rhs);
        rhs_operators[relative_level].initialize(solve_group,
                                                 solutions.get_matrix_free(
                                                   relative_level),
                                                 solutions.get_global_to_block_index());
      }
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
                            attributes_list[global_index].name);
                if (name_it != initial_condition_file.simulation_variable_names.end())
                  {
                    dealii::VectorTools::interpolate(
                      solve_context->get_mapping(),
                      solve_context->get_dof_manager().get_dof_handler(global_index),
                      ReadInitialCondition<dim, number>(
                        *name_it,
                        attributes_list[global_index].field_type,
                        initial_condition_file,
                        solve_context->get_user_inputs().get_spatial_discretization()),
                      solutions.get_solution_vector(global_index));
                  }
              }
          }
        else
          {
            dealii::VectorTools::interpolate(
              solve_context->get_mapping(),
              solve_context->get_dof_manager().get_dof_handler(index),
              InitialCondition<dim, degree, number>(
                index,
                attributes_list[global_index].field_type,
                solve_context->get_pde_operator()),
              solutions.get_solution_vector(global_index));
          }
        solutions.apply_initial_condition_for_old_fields();
      }
  }

  /**
   * @brief Get the invm handler.
   */
  [[nodiscard]] const InvmHandler<dim, degree, number> &
  get_invm_handler() const
  {
    return solve_context->get_invm_handler();
  }

  /**
   * @brief Get the solution handler.
   */
  [[nodiscard]] GroupSolutionHandler<dim, number> &
  get_solution_manager() const;

  /**
   * @brief Get the element volume container.
   */
  [[nodiscard]] const ElementVolumeContainer<dim, degree, number> &
  get_element_volume_container() const
  {
    return solve_context->get_element_volume_container();
  }

  /**
   * @brief Get the pde operator.
   */
  [[nodiscard]] const std::shared_ptr<const PDEOperator<dim, degree, number>> &
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
   * @brief List of field attributes available.
   */
  const std::vector<FieldAttributes> attributes_list;
  /**
   * @brief Information about the solve group this handler is responsible for.
   */
  SolveGroup solve_group;
  /**
   * @brief Solution vectors for fields handled by this solver.
   */
  GroupSolutionHandler<dim, number> solutions;

  std::vector<MFOperator<dim, degree, number>> rhs_operators;

  std::vector<GroupSolverBase<dim, degree, number> *> aux_solvers;

  /**
   * @brief Solver context provides access to external information.
   */
  const std::shared_ptr<SolveContext<dim, degree, number>> solve_context;
};

PRISMS_PF_END_NAMESPACE
