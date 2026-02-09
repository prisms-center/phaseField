// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/types.h>

#include <prismspf/solvers/solvers.h>

#include <prismspf/utilities/integrator.h>

#include <prismspf/config.h>
#include <prismspf/nucleation/nucleation.h>
#include <prismspf/nucleation/nucleus_refinement_function.h>

#include "prismspf/core/field_attributes.h"
#include "prismspf/core/group_solution_handler.h"

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
  Problem(const std::vector<FieldAttributes>                            &field_attributes,
          const std::vector<SolveGroup>                                 &solve_groups,
          const UserInputParameters<dim>                                &_user_inputs,
          PhaseFieldTools<dim>                                          &_pf_tools,
          const std::shared_ptr<const PDEOperator<dim, degree, number>> &_pde_operator);

  /**
   * @brief Run initialization and solving steps of the given problem.
   */
  void
  run();

  /**
   * @brief Scalar FESystem.
   */
  static const dealii::FESystem<dim>
    scalar_fe_system(dealii::FE_Q<dim>(dealii::QGaussLobatto<1>(degree + 1)), 1);
  /**
   * @brief Scalar FESystem.
   */
  static const dealii::FESystem<dim>
    vector_fe_system(dealii::FE_Q<dim>(dealii::QGaussLobatto<1>(degree + 1)), dim);
  static const std::array<const dealii::FESystem<dim>, 2> fe_systems;

private:
  /**
   * @brief Main time-stepping loop that calls solve_increment, reinit_system,
   * output_results, etc...
   */
  void
  solve();

  /**
   * @brief Solve a single increment of the given PDEs.
   */
  void
  solve_increment();

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
   * @brief Solvers.
   */
  std::vector<std::shared_ptr<GroupSolverBase<dim, degree, number>>> solvers;

  /**
   * @brief Triangulation handler.
   */
  TriangulationManager<dim> triangulation_handler;

  /**
   * @brief Constraint handler.
   */
  ConstraintHandler<dim, degree, number> constraint_handler;

  /**
   * @brief Solution handler.
   */
  GroupSolutionHandler<dim, number> solution_handler;

  /**
   * @brief DoF manager.
   */
  DofManager<dim> dof_manager;

  /**
   * @brief Collection of finite element systems. This is just a collection of two
   * FESystem's: one for scalar fields and one for vector fields. For now they both use
   * FE_Q finite elements.
   */
  std::map<FieldInfo::TensorRank, dealii::FESystem<dim>> fe_system;

  /**
   * @brief Mappings to and from reference cell.
   */
  dealii::MappingQ1<dim> mapping;

  /**
   * @brief Solver context.
   */
  SolverContext<dim, degree, number> solver_context;

  /**
   * @brief Grid refiner.
   */
  GridRefiner<dim, degree, number> grid_refiner;
};

PRISMS_PF_END_NAMESPACE
