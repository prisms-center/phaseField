// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/refinement_manager.h>
#include <prismspf/core/simulation_timer.h>
#include <prismspf/core/types.h>

#include <prismspf/nucleation/nucleation_manager.h>
#include <prismspf/nucleation/nucleus_refinement_function.h>

#include <prismspf/solvers/solve_context.h>

#include <prismspf/utilities/integrator.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This is the main class that handles the construction and solving of
 * user-specified PDEs.
 */
template <unsigned int dim, unsigned int degree, typename number>
class Problem
{
public:
  /**
   * @brief Constructor.
   */
  Problem(const std::vector<FieldAttributes>                          &field_attributes,
          const std::vector<SolveGroup>                               &solve_groups,
          const UserInputParameters<dim>                              &_user_inputs,
          PhaseFieldTools<dim>                                        &_pf_tools,
          const std::shared_ptr<PDEOperatorBase<dim, degree, number>> &_pde_operator);

  /**
   * @brief Main time-stepping loop that calls solve_increment, reinit_system,
   * output_results, etc...
   */
  void
  solve();

private:
  /**
   * @brief Solve a single increment of the given PDEs.
   */
  void
  solve_increment(SimulationTimer &sim_timer);

  /**
   * @brief Initialize the system.
   */
  void
  init_system();

  /**
   * @brief Reinitialize the system.
   */
  void
  reinit_system();

  /**
   * @brief Field attributes.
   */
  std::vector<FieldAttributes> field_attributes;

  /**
   * @brief Solve groups.
   */
  std::vector<SolveGroup> solve_groups;

  /**
   * @brief User-inputs.
   */
  const UserInputParameters<dim> *user_inputs_ptr;

  /**
   * @brief Phase field tools.
   */
  PhaseFieldTools<dim> *pf_tools;

  /**
   * @brief Triangulation handler.
   */
  TriangulationManager<dim> triangulation_manager;

  /**
   * @brief DoF manager.
   */
  DofManager<dim> dof_manager;

  /**
   * @brief Constraint handler.
   */
  ConstraintManager<dim, degree, number> constraint_manager;

  /**
   * @brief Solvers.
   */
  std::vector<std::shared_ptr<GroupSolverBase<dim, degree, number>>> solvers;

  /**
   * @brief Solution indexer
   */
  SolutionIndexer<dim, number> solution_indexer;

  /**
   * @brief Solver context.
   */
  SolveContext<dim, degree, number> solve_context;

  /**
   * @brief Grid refiner.
   */
  RefinementManager<dim, degree, number> grid_refiner;
};

PRISMS_PF_END_NAMESPACE
