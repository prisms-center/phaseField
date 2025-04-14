// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/pde_operator.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

#ifdef PRISMS_PF_WITH_CALIPER
#  include <caliper/cali.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
class userInputParameters;

template <int dim>
class constraintHandler;

template <int dim>
class dofHandler;

template <int dim, int degree, typename number>
class invmHandler;

template <int dim, typename number>
class matrixfreeHandler;

template <int dim>
class solutionHandler;

template <int dim, typename number>
class solutionOutput;

class timer;

template <int dim>
class triangulationHandler;

struct variableAttributes;

template <int dim, int degree>
class explicitConstantSolver;

template <int dim, int degree>
class explicitPostprocessSolver;

template <int dim, int degree>
class explicitSolver;

template <int dim, int degree>
class nonexplicitAuxiliarySolver;

template <int dim, int degree>
class nonexplicitLinearSolver;

template <int dim, int degree>
class nonexplicitSelfNonlinearSolver;

template <int dim, int degree, typename number>
class computeIntegral;

template <int dim, int degree, typename number>
class elementVolume;

/**
 * \brief This is the main class that handles the construction and solving of
 * user-specified PDEs.
 */
template <int dim, int degree>
class PDEProblem
{
public:
  /**
   * \brief Constructor.
   */
  explicit PDEProblem(
    const userInputParameters<dim>                         &_user_inputs,
    std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator,
    std::shared_ptr<const PDEOperator<dim, degree, float>>  _pde_operator_float);

  /**
   * \brief Run initialization and solving steps of the given problem.
   */
  void
  run();

private:
  /**
   * \brief Main time-stepping loop that calls solve_increment, reinit_system,
   * output_results, etc...
   */
  void
  solve();

  /**
   * \brief Solve a single increment of the given PDEs.
   */
  void
  solve_increment();

  /**
   * \brief Initialize the system.
   */
  void
  init_system();

  /**
   * \brief Reinitialize the system.
   */
  void
  reinit_system();

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> *user_inputs;

  /**
   * \brief Triangulation handler.
   */
  triangulationHandler<dim> triangulation_handler;

  /**
   * \brief Constraint handler.
   */
  constraintHandler<dim> constraint_handler;

  /**
   * \brief Matrix-free object handler for non-multigrid data.
   */
  matrixfreeHandler<dim, double> matrix_free_handler;

  /**
   * \brief Matrix-free object handler for multigrid data.
   */
  dealii::MGLevelObject<matrixfreeHandler<dim, float>> multigrid_matrix_free_handler;

  /**
   * \brief invm handler.
   */
  invmHandler<dim, degree, double> invm_handler;

  /**
   * \brief Solution handler.
   */
  solutionHandler<dim> solution_handler;

  /**
   * \brief DoF handler.
   */
  dofHandler<dim> dof_handler;

  /**
   * \brief Collection of finite element systems. This is just a collection of two
   * FESystem's: one for scalar fields and one for vector fields. For now they both use
   * FE_Q finite elements.
   */
  std::map<fieldType, dealii::FESystem<dim>> fe_system;

  /**
   * \brief Mappings to and from reference cell.
   */
  dealii::MappingQ1<dim> mapping;

  /**
   * \brief Element volumes.
   */
  elementVolume<dim, degree, double> element_volume;

  /**
   * \brief Integral utility.
   *
   * TODO (landinjm): Rename this class.
   */
  computeIntegral<dim, degree, double> integral_computer;

  /**
   * \brief Explicit constant field solver class.
   */
  explicitConstantSolver<dim, degree> explicit_constant_solver;

  /**
   * \brief Explicit field solver class.
   */
  explicitSolver<dim, degree> explicit_solver;

  /**
   * \brief Postprocessed explicit field solver class.
   */
  explicitPostprocessSolver<dim, degree> postprocess_explicit_solver;

  /**
   * \brief Nonexplicit auxiliary field solver class.
   */
  nonexplicitAuxiliarySolver<dim, degree> nonexplicit_auxiliary_solver;

  /**
   * \brief Nonexplicit linear field solver class.
   */
  nonexplicitLinearSolver<dim, degree> nonexplicit_linear_solver;

  /**
   * \brief Nonexplicit self nonlinear field solver class.
   */
  nonexplicitSelfNonlinearSolver<dim, degree> nonexplicit_self_nonlinear_solver;
};

PRISMS_PF_END_NAMESPACE
