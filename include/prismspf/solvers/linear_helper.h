// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/invm_manager.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/mf_operator.h>
#include <prismspf/solvers/solver_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class exists as a hack to get access to the residual vector for the custom
 * convergence criterion because dealii doesn't provide it for some reason. The function
 * `print_vectors` gets called before the convergence check in the same scope and is
 * overridable.
 * Be careful, because the pointers are only valid during the linear solve, so they
 * shouldn't be called outside of that context.
 */
template <typename number>
class CGSolver : public dealii::SolverCG<BlockVector<number>>
{
public:
  using AddData   = dealii::SolverCG<BlockVector<number>>::AdditionalData;
  using BaseClass = dealii::SolverCG<BlockVector<number>>;
  using BaseClass::all_condition_numbers_signal;
  using BaseClass::all_eigenvalues_signal;
  using BaseClass::condition_number_signal;
  using BaseClass::determine_beta_by_flexible_formula;
  using BaseClass::eigenvalues_signal;

  /**
   * @brief Constructor.
   */
  CGSolver(dealii::SolverControl                     &cn,
           dealii::VectorMemory<BlockVector<number>> &mem,
           const AddData                             &data = AddData())
    : dealii::SolverCG<BlockVector<number>>(cn, mem, data)
  {}

  /**
   * @brief Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  explicit CGSolver(dealii::SolverControl &cn, const AddData &data = AddData())
    : dealii::SolverCG<BlockVector<number>>(cn, data)
  {}

  /**
   * @brief Solve method we are replacing in prisms-pf because we need access to the
   * residual.
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  solve_better(const MatrixType          &A,
               BlockVector<number>       &x,
               const BlockVector<number> &b,
               const PreconditionerType  &preconditioner)
  {
    using namespace dealii;
    using VectorType = BlockVector<number>;

    SolverControl::State solver_state = SolverControl::iterate;

    LogStream::Prefix prefix("cg");

    // Should we build the matrix for eigenvalue computations?
    const bool do_eigenvalues =
      !condition_number_signal.empty() || !all_condition_numbers_signal.empty() ||
      !eigenvalues_signal.empty() || !all_eigenvalues_signal.empty();

    // vectors used for eigenvalue computations
    std::vector<typename VectorType::value_type> diagonal;
    std::vector<typename VectorType::value_type> offdiagonal;

    typename VectorType::value_type eigen_beta_alpha = 0;

    int it = 0;

    internal::SolverCG::IterationWorker<VectorType, MatrixType, PreconditionerType>
      worker(A, preconditioner, determine_beta_by_flexible_formula, this->memory, x);

    worker.startup(b);

    solver_state = this->iteration_status(0, worker.residual_norm, x);
    if (solver_state != SolverControl::iterate)
      {
        return;
      }

    while (solver_state == SolverControl::iterate)
      {
        ++it;

        worker.do_iteration(it);

        print_vectors(it, x, worker.r, worker.p);

        if (it > 1)
          {
            this->coefficients_signal(worker.previous_alpha, worker.beta);
            // set up the vectors containing the diagonal and the off diagonal
            // of the projected matrix.
            if (do_eigenvalues)
              {
                diagonal.push_back(number(1.) / worker.previous_alpha + eigen_beta_alpha);
                eigen_beta_alpha = worker.beta / worker.previous_alpha;
                offdiagonal.push_back(std::sqrt(worker.beta) / worker.previous_alpha);
              }
            compute_eigs_and_cond(diagonal,
                                  offdiagonal,
                                  all_eigenvalues_signal,
                                  all_condition_numbers_signal);
          }

        solver_state = this->iteration_status(it, worker.residual_norm, x);
      }

    worker.finalize_after_convergence(it);

    compute_eigs_and_cond(diagonal,
                          offdiagonal,
                          eigenvalues_signal,
                          condition_number_signal);

    AssertThrow(solver_state == SolverControl::success,
                SolverControl::NoConvergence(it, worker.residual_norm));
  }
};

// in linearsolver
/* dealii::SolverControl::State
rmse(const BlockVector<number>              &residual,
     number                                  tolerance,
     const InvMManager<dim, degree, number> &invm,
     unsigned int                            relative_level)
{
  const auto &field_to_block = solutions.get_global_to_block_index();
  number      sum            = 0.0;
  for (unsigned int field_index : solve_block.field_indices)
    {
      const unsigned int block_index = field_to_block[field_index];
      TensorRank rank = solve_context->get_field_attributes()[field_index].field_type;

      SolutionVector<number> vec = residual.block(block_index); // residual
      vec.scale(vec);                                           // l2

      vec.scale(invm.get_jxw(rank, relative_level));
      sum += vec.mean_value();
    }
  using std::sqrt;
  if (sqrt(sum) < tolerance)
    {
      return dealii::SolverControl::success;
    }
  return dealii::SolverControl::iterate;
} */

// before solving
/* cg_solver.connect(
  [&](const unsigned int, const double, const BlockVector<number> &)
  {
    return rmse(*(cg_solver.r_ptr),
                params.tolerance,
                solve_context->get_invm_manager(),
                relative_level);
  }); */
