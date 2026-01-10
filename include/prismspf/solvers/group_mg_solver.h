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

#include <prismspf/core/timer.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/group_linear_solver.h>
#include <prismspf/solvers/mf_operator.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolverContext;

/**
 * @brief This class handles the explicit solves of all explicit fields
 */
template <unsigned int dim, unsigned int degree, typename number>
class MGSolver : public LinearSolver<dim, degree, number>
{
  using GroupSolverBase = GroupSolverBase<dim, degree, number>;
  using LinearSolver    = LinearSolver<dim, degree, number>;
  using GroupSolverBase::rhs_operators;
  using GroupSolverBase::solutions;
  using GroupSolverBase::solve_group;
  using GroupSolverBase::solver_context;
  using LinearSolver::do_linear_solve;
  using LinearSolver::lhs_operators;
  using BlockVector = SolutionHandler<dim, number>::BlockVector;

public:
  /**
   * @brief Constructor.
   */
  MGSolver(SolveGroup                                _solve_group,
           const SolverContext<dim, degree, number> &_solver_context)
    : LinearSolver(_solve_group, _solver_context)
  {}

  /**
   * @brief Solve for a single update step.
   */
  void
  solve_level(unsigned int relative_level) override
  {}

  void
  mg_solve()
  {
    int min_level = 0;
    int max_level = 0;
    // 1. Level operators
    dealii::MGLevelObject<MFOperator<dim, degree, number>> mg_lhs_operators(min_level,
                                                                            max_level);
    // TODO: Initialize them

    // 2. Solution and rhs vectors (rhs and solution here are used as scratch space by mg,
    // but are not actually the same as the top level vectors)
    dealii::MGLevelObject<BlockVector> mg_solution_vectors(min_level, max_level);
    dealii::MGLevelObject<BlockVector> mg_rhs_vectors(min_level, max_level);

    // TODO: Initialize their shapes

    // 3. MG transfer
    dealii::MGTransferBlockMatrixFree<dim, number> mg_transfer; // Constraints?
    // NOTE: dof_handler.distribute_mg_dofs() must have been called
    mg_transfer.build(
      solver_context->dof_manager.get_dof_handlers(solve_group.field_indices, 0));

    // 4. MG Smoother (takes in operators) This is similar to a solver, but is
    // conceptually different.
    using SmootherPrecond =
      dealii::PreconditionChebyshev<MFOperator<dim, degree, number>, BlockVector>;
    dealii::MGLevelObject<typename SmootherPrecond::AdditionalData> smoother_data(
      min_level,
      max_level);
    for (unsigned int level = min_level; level <= max_level; ++level)
      {
        // smoother_data[level].smoothing_range;
        // smoother_data[level].degree;
        // smoother_data[level].eig_cg_n_iterations;
        // smoother_data[level].preconditioner;
        // smoother_data[level].constraints;
      }
    // Wrapper around a generic preconditioner to be used as a smoother
    using Smoother = dealii::MGSmootherPrecondition<MFOperator<dim, degree, number>,
                                                    SmootherPrecond,
                                                    BlockVector>;
    Smoother mg_smoother;
    mg_smoother.initialize(mg_lhs_operators, smoother_data);

    // The MG method requires the choice of a solve on the coarsest level. There are a
    // couple of main options:
    // A. Direct solve // we're not doing this
    // B. Iterative solve (using a Krylov method, (like CG))
    /* dealii::MGCoarseGridIterativeSolver<BlockVector,
                                        dealii::SolverCG<BlockVector>,
                                        MFOperator<dim, degree, number>,
                                        dealii::PreconditionIdentity> mg_coarse_solver;*/
    // C. Smoothing operation (don't actually solve, just use the same smoother as on
    // other levels)
    /* dealii::MGCoarseGridApplySmoother<BlockVector> mg_coarse_solver;
    mg_coarse_solver.initialize(mg_smoother); */

    // 5. Coarse grid solver
    dealii::MGCoarseGridApplySmoother<BlockVector> mg_coarse_solver;
    mg_coarse_solver.initialize(mg_smoother);

    // 6. Multigrid object
    dealii::Multigrid<BlockVector> multigrid(
      mg_lhs_operators,
      mg_coarse_solver,
      mg_transfer,
      mg_smoother,
      mg_smoother,
      min_level,
      max_level,
      dealii::Multigrid<BlockVector>::Cycle::v_cycle);

    // 7. Turn MG into a preconditioner object
    dealii::
      PreconditionMG<dim, BlockVector, dealii::MGTransferBlockMatrixFree<dim, number>>
        preconditioner_mg(
          solver_context->dof_manager.get_dof_handlers(solve_group.field_indices, 0),
          multigrid,
          mg_transfer);

    // 8. Solve the system with CG + MG preconditioner
    dealii::SolverControl         cg_solver_control(/* TODO */);
    dealii::SolverCG<BlockVector> cg_solver(cg_solver_control);
  }

private:
  /**
   * @brief Solver control. Contains max iterations and tolerance.
   */
  dealii::SolverControl linear_solver_control;
};

PRISMS_PF_END_NAMESPACE
