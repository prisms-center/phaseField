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

#include <prismspf/core/multigrid_info.h>

#include <prismspf/solvers/linear_solver_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class DofHandler;

/**
 * @brief Class that handles the assembly and solving of a field with a GMG preconditioner
 */
template <unsigned int dim, unsigned int degree, typename number>
class GMGSolver : public LinearSolverBase<dim, degree, number>
{
public:
  using SystemMatrixType = MatrixFreeOperator<dim, degree, number>;
  using LevelMatrixType  = MatrixFreeOperator<dim, degree, float>;
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<number>;
  using MGVectorType     = dealii::LinearAlgebra::distributed::Vector<float>;

  /**
   * @brief Constructor.
   */
  GMGSolver(const UserInputParameters<dim>                       &_user_inputs,
            const VariableAttributes                             &_variable_attributes,
            const MatrixfreeHandler<dim, number>                 &_matrix_free_handler,
            const ConstraintHandler<dim, degree, number>         &_constraint_handler,
            const TriangulationHandler<dim>                      &_triangulation_handler,
            const DofHandler<dim>                                &_dof_handler,
            dealii::MGLevelObject<MatrixfreeHandler<dim, float>> &_mg_matrix_free_handler,
            SolutionHandler<dim, number>                         &_solution_handler,
            std::shared_ptr<const PDEOperator<dim, degree, number>> _pde_operator,
            std::shared_ptr<const PDEOperator<dim, degree, float>>  _pde_operator_float,
            const MGInfo<dim>                                      &_mg_info);

  /**
   * @brief Destructor.
   */
  ~GMGSolver() override = default;

  /**
   * @brief Initialize the system.
   */
  void
  init() override;

  /**
   * @brief Reinitialize the system.
   */
  void
  reinit() override;

  /**
   * @brief Solve the system Ax=b.
   */
  void
  solve(const number &step_length = 1.0) override;

private:
  /**
   * @brief Triangulation handler.
   */
  const TriangulationHandler<dim> *triangulation_handler;

  /**
   * @brief DoF handler.
   */
  const DofHandler<dim> *dof_handler;

  /**
   * @brief Matrix-free object handler for multigrid data.
   */
  dealii::MGLevelObject<MatrixfreeHandler<dim, float>> *mg_matrix_free_handler;

  /**
   * @brief PDE operator but for floats!
   */
  std::shared_ptr<const PDEOperator<dim, degree, float>> pde_operator_float;

  /**
   * @brief Multigrid information.
   */
  const MGInfo<dim> *mg_info;

  /**
   * @brief Minimum multigrid level
   */
  unsigned int min_level = 0;

  /**
   * @brief Maximum multigrid level
   */
  unsigned int max_level = 0;

  /**
   * @brief The local index of the change variable.
   */
  Types::Index change_local_index = Numbers::invalid_index;

  /**
   * @brief Mappings to and from reference cell.
   *
   * TODO (landinjm): This should be the same as the rest of the problem.
   */
  dealii::MappingQ1<dim> mapping;

  /**
   * @brief Collection of transfer operators for each multigrid level.
   */
  std::vector<dealii::MGLevelObject<dealii::MGTwoLevelTransfer<dim, MGVectorType>>>
    mg_transfer_operators;

  /**
   * @brief PDE operator for each multigrid level.
   */
  std::unique_ptr<dealii::MGLevelObject<LevelMatrixType>> mg_operators;

  /**
   * @brief Multigrid object for storing all operators.
   */
  std::shared_ptr<dealii::mg::Matrix<MGVectorType>> mg_matrix;

  /**
   * @brief Transfer operator for global coarsening.
   */
  std::vector<std::shared_ptr<dealii::MGTransferGlobalCoarsening<dim, MGVectorType>>>
    mg_transfer;
};

PRISMS_PF_END_NAMESPACE
