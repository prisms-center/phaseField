// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <prismspf/core/dof_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_output.h>
#include <prismspf/core/triangulation_handler.h>

#include <prismspf/solvers/linear_solver_base.h>

#include <prismspf/config.h>

#ifdef PRISMS_PF_WITH_CALIPER
#  include <caliper/cali.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Class that handles the assembly and solving of a field with a GMG preconditioner
 */
template <int dim, int degree>
class GMGSolver : public linearSolverBase<dim, degree>
{
public:
  using SystemMatrixType = matrixFreeOperator<dim, degree, double>;
  using LevelMatrixType  = matrixFreeOperator<dim, degree, float>;
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<double>;
  using MGVectorType     = dealii::LinearAlgebra::distributed::Vector<float>;

  /**
   * \brief Constructor.
   */
  GMGSolver(const userInputParameters<dim>                       &_user_inputs,
            const variableAttributes                             &_variable_attributes,
            const matrixfreeHandler<dim>                         &_matrix_free_handler,
            const constraintHandler<dim>                         &_constraint_handler,
            const triangulationHandler<dim>                      &_triangulation_handler,
            const dofHandler<dim>                                &_dof_handler,
            dealii::MGLevelObject<matrixfreeHandler<dim, float>> &_mg_matrix_free_handler,
            solutionHandler<dim>                                 &_solution_handler,
            std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator,
            std::shared_ptr<const PDEOperator<dim, degree, float>>  _pde_operator_float);

  /**
   * \brief Destructor.
   */
  ~GMGSolver() override;

  /**
   * \brief Initialize the system.
   */
  void
  init() override;

  /**
   * \brief Reinitialize the system.
   */
  void
  reinit() override;

  /**
   * \brief Solve the system Ax=b.
   */
  void
  solve(const double &step_length = 1.0) override;

private:
  /**
   * \brief Triangulation handler.
   */
  const triangulationHandler<dim> *triangulation_handler;

  /**
   * \brief DoF handler.
   */
  const dofHandler<dim> *dof_handler;

  /**
   * \brief Matrix-free object handler for multigrid data.
   */
  dealii::MGLevelObject<matrixfreeHandler<dim, float>> *mg_matrix_free_handler;

  /**
   * \brief Minimum multigrid level
   */
  unsigned int min_level = 0;

  /**
   * \brief Maximum multigrid level
   */
  unsigned int max_level = 0;

  /**
   * \brief Mappings to and from reference cell.
   *
   * TODO (landinjm): This should be the same as the rest of the problem.
   */
  dealii::MappingQ1<dim> mapping;

  /**
   * \brief Collection of transfer operators for each multigrid level.
   */
  dealii::MGLevelObject<dealii::MGTwoLevelTransfer<dim, MGVectorType>>
    mg_transfer_operators;

  /**
   * \brief PDE operator for each multigrid level.
   */
  std::unique_ptr<dealii::MGLevelObject<LevelMatrixType>> mg_operators;

  /**
   * \brief Multigrid object for storing all operators.
   */
  std::shared_ptr<dealii::mg::Matrix<MGVectorType>> mg_matrix;

  /**
   * \brief Transfer operator for global coarsening.
   */
  std::shared_ptr<dealii::MGTransferGlobalCoarsening<dim, MGVectorType>> mg_transfer;

  /**
   * \brief Subset of fields that are necessary for the source of the newton update for
   * each multigrid level.
   */
  dealii::MGLevelObject<std::vector<MGVectorType *>> mg_newton_update_src;

  /**
   * \brief PDE operator but for floats!
   */
  std::shared_ptr<const PDEOperator<dim, degree, float>> pde_operator_float;
};

PRISMS_PF_END_NAMESPACE