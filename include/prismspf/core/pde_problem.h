// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/grid_refiner.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/solver_handler.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/element_volume.h>
#include <prismspf/utilities/integrator.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This is the main class that handles the construction and solving of
 * user-specified PDEs.
 */
template <unsigned int dim, unsigned int degree, typename number>
class PDEProblem
{
public:
  /**
   * @brief Constructor.
   */
  PDEProblem(
    const UserInputParameters<dim>                                &_user_inputs,
    const std::shared_ptr<const PDEOperator<dim, degree, number>> &_pde_operator,
    const std::shared_ptr<const PDEOperator<dim, degree, float>>  &_pde_operator_float);

  /**
   * @brief Run initialization and solving steps of the given problem.
   */
  void
  run();

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
   * @brief User-inputs.
   */
  const UserInputParameters<dim> *user_inputs;

  /**
   * @brief Multigrid info class.
   */
  MGInfo<dim> mg_info;

  /**
   * @brief Triangulation handler.
   */
  TriangulationHandler<dim> triangulation_handler;

  /**
   * @brief Constraint handler.
   */
  ConstraintHandler<dim, degree, number> constraint_handler;

  /**
   * @brief Matrix-free object handler for non-multigrid data.
   */
  MatrixfreeHandler<dim, number> matrix_free_handler;

  /**
   * @brief Matrix-free object handler for multigrid data.
   */
  dealii::MGLevelObject<MatrixfreeHandler<dim, float>> multigrid_matrix_free_handler;

  /**
   * @brief invm handler.
   */
  InvmHandler<dim, degree, number> invm_handler;

  /**
   * @brief Solution handler.
   */
  SolutionHandler<dim, number> solution_handler;

  /**
   * @brief DoF handler.
   */
  DofHandler<dim> dof_handler;

  /**
   * @brief Collection of finite element systems. This is just a collection of two
   * FESystem's: one for scalar fields and one for vector fields. For now they both use
   * FE_Q finite elements.
   */
  std::map<FieldType, dealii::FESystem<dim>> fe_system;

  /**
   * @brief Mappings to and from reference cell.
   */
  dealii::MappingQ1<dim> mapping;

  /**
   * @brief Solver context.
   */
  SolverContext<dim, degree, number> solver_context;

  /**
   * @brief Element volumes.
   */
  ElementVolume<dim, degree, number> element_volume;

  /**
   * @brief Integrator utility.
   */
  Integrator<dim, degree, number> integrator;

  /**
   * @brief Grid refiner context.
   */
  GridRefinementContext<dim, degree, number> grid_refiner_context;

  /**
   * @brief Grid refiner.
   */
  GridRefiner<dim, degree, number> grid_refiner;

  /**
   * @brief Solver handler.
   */
  SolverHandler<dim, degree, number> solver_handler;
};

PRISMS_PF_END_NAMESPACE
