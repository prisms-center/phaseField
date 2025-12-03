// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/types.h>
#include <prismspf/core/vector_mapping.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

template <unsigned int dim>
struct PhaseFieldTools;

template <unsigned int dim, unsigned int degree, typename number>
class PDEOperator;

template <unsigned int dim, unsigned int degree, typename number>
class ConstraintHandler;

template <unsigned int dim>
class DofHandler;

template <unsigned int dim, unsigned int degree, typename number>
class GridRefiner;

template <unsigned int dim, unsigned int degree, typename number>
class InvmHandler;

template <unsigned int dim, typename number>
class MatrixFreeContainer;

template <unsigned int dim>
class MGInfo;

template <unsigned int dim, typename number>
class SolutionHandler;

template <unsigned int dim, unsigned int degree, typename number>
class SolverHandler;

template <unsigned int dim>
class TriangulationHandler;

template <unsigned int dim, unsigned int degree, typename number>
class GridRefinementContext;

template <unsigned int dim, unsigned int degree, typename number>
class SolverContext;

template <unsigned int dim, unsigned int degree, typename number>
class ElementVolumeContainer;

template <unsigned int dim, unsigned int degree, typename number>
class Integrator;

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
    PhaseFieldTools<dim>                                          &_pf_tools,
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
   * @brief Phase field tools.
   */
  PhaseFieldTools<dim> *pf_tools;

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
   * @brief Matrix-free object container.
   */
  MatrixFreeContainer<dim, number> matrix_free_container;

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
  std::map<FieldInfo::TensorRank, dealii::FESystem<dim>> fe_system;

  /**
   * @brief Mappings to and from reference cell.
   */
  dealii::MappingQ1<dim> mapping;

  /**
   * @brief Element volume container.
   */
  ElementVolumeContainer<dim, degree, number> element_volume_container;

  /**
   * @brief Solver context.
   */
  SolverContext<dim, degree, number> solver_context;

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
