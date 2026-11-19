// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_selector.h>

#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/linear_solver.h>
#include <prismspf/solvers/mf_operator.h>
#include <prismspf/solvers/solver_base.h>

#include <prismspf/user_inputs/solve_parameters.h>

#include <prismspf/utilities/timer.h>

#include <prismspf/config.h>

#include <cmath>
#include <memory>
#include <unordered_set>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolveContext;

/**
 * @brief Matrix-free operator wrapper that enforces homogeneous active-set
 * constraints (Delta u_i = 0 for i in active_set) on the Jacobian operator
 * without rebuilding MatrixFree every Newton step.
 */
template <unsigned int dim, unsigned int degree, typename number>
class ActiveSetOperatorWrapper : public MATRIX_FREE_OPERATOR_BASE
{
public:
  ActiveSetOperatorWrapper(MFOperator<dim, degree, number> &base_op,
                           const unsigned int               constrained_block,
                           const dealii::IndexSet          &active_set)
    : base_operator(base_op)
    , block_idx(constrained_block)
    , active_dofs(active_set)
  {}

  dealii::types::global_dof_index
  m() const
  {
    return base_operator.m();
  }

  number
  el(const unsigned int &row, const unsigned int &col) const
  {
    return base_operator.el(row, col);
  }

  void
  vmult(BlockVector<number> &dst, const BlockVector<number> &src) const
  {
    BlockVector<number> filtered_src(src);
    filtered_src.zero_out_ghost_values();
    for (const auto idx : active_dofs)
      {
        if (filtered_src.block(block_idx).in_local_range(idx))
          {
            filtered_src.block(block_idx)[idx] = 0.0;
          }
      }
    filtered_src.update_ghost_values();

    base_operator.vmult(dst, filtered_src);

    // Pin active rows to identity (dst_i = src_i) and keep dst non-ghosted
    // as required by MFOperator::vmult.
    for (const auto idx : active_dofs)
      {
        if (dst.block(block_idx).in_local_range(idx))
          {
            dst.block(block_idx)[idx] = src.block(block_idx)[idx];
          }
      }
  }

private:
  MFOperator<dim, degree, number> &base_operator;
  const unsigned int               block_idx;
  const dealii::IndexSet          &active_dofs;
};

/**
 * @brief Nonlinear Newton solver with a primal-dual active-set strategy
 * (modeled after deal.II step-42) for obstacle / variational inequality problems.
 */
template <unsigned int dim, unsigned int degree, typename number>
class ActiveSetNewtonSolver : public LinearSolver<dim, degree, number>
{
protected:
  using SolverBase<dim, degree, number>::solutions;
  using SolverBase<dim, degree, number>::solve_context;
  using SolverBase<dim, degree, number>::solve_block;
  using LinearSolver<dim, degree, number>::normalization_value;
  using LinearSolver<dim, degree, number>::lhs_operator;
  using LinearSolver<dim, degree, number>::rhs_operator;
  using LinearSolver<dim, degree, number>::rhs_vector;

  using ActiveSetPreconditionChebyshev =
    dealii::PreconditionChebyshev<ActiveSetOperatorWrapper<dim, degree, number>,
                                  BlockVector<number>,
                                  dealii::DiagonalMatrix<BlockVector<number>>>;

public:
  ActiveSetNewtonSolver(SolveBlock                               _solve_block,
                        const SolveContext<dim, degree, number> &_solve_context,
                        const unsigned int _constrained_block = 0,
                        const number       _penalty_c         = 1.0e4)
    : LinearSolver<dim, degree, number>(_solve_block, _solve_context)
    , constrained_block(_constrained_block)
    , penalty_c(_penalty_c)
  {}

  void
  init(const std::list<SolveBlock> &all_solve_blocks) override
  {
    LinearSolver<dim, degree, number>::init(all_solve_blocks);
    initialize_active_set_linear_solver();
    setup_active_set_structures();
  }

  void
  reinit() override
  {
    LinearSolver<dim, degree, number>::reinit();
    initialize_active_set_linear_solver();
    setup_active_set_structures();
  }

  void
  solve_impl() override
  {
    Timer::Scope scope("Active-Set Nonlinear Solver");

    const number       base_step_length      = newton_params().step_length;
    const number       newton_tolerance      = newton_params().tolerance_value * normalization_value();
    const unsigned int newton_max_iterations = newton_params().max_iterations;

    BlockVector<number>             &newton_residual = rhs_vector;
    MFOperator<dim, degree, number> &rhs_op          = rhs_operator;
    MFOperator<dim, degree, number> &lhs_op          = lhs_operator;

    ActiveSetOperatorWrapper<dim, degree, number> constrained_lhs_op(
      lhs_op, constrained_block, active_set);

    if (!solutions.get_primary_solutions().old_solutions.empty())
      {
        solutions.get_solution_full_vector() = solutions.get_old_solution_full_vector(0);
      }

    active_set.clear();
    dealii::IndexSet old_active_set(active_set);

    bool         newton_unconverged = true;
    unsigned int iter               = 0;
    number       l2_norm            = -1.0;
    int          total_lin_iters    = 0;

    while (newton_unconverged && iter < newton_max_iterations)
      {
        // 1. Enforce base Dirichlet and hanging-node constraints
        solutions.apply_constraints();
        solutions.update_ghosts();

        // 2. Compute UNCONDENSED residual (without active-set zeroing)
        //    so reaction forces / multipliers lambda_i are preserved.
        newton_residual.zero_out_ghost_values();
        rhs_op.compute_operator(newton_residual);
        newton_residual.update_ghost_values();

        // 3. Update active set and project active DoFs onto the obstacle
        old_active_set = active_set;
        update_active_set_and_project(newton_residual);

        // Re-evaluate residual at the newly projected state if active set changed
        const unsigned int active_set_diff = dealii::Utilities::MPI::sum(
          (active_set == old_active_set) ? 0u : 1u, MPI_COMM_WORLD);

        if (active_set_diff > 0)
          {
            solutions.apply_constraints();
            solutions.update_ghosts();

            newton_residual.zero_out_ghost_values();
            rhs_op.compute_operator(newton_residual);
            newton_residual.update_ghost_values();
          }

        // 4. Condense active DoFs out of the residual for the linear solve
        //    and convergence check (residual on active DoFs = reaction force).
        zero_active_dofs(newton_residual);

        l2_norm = newton_residual.l2_norm();
        newton_unconverged =
          (l2_norm > newton_tolerance) || (active_set_diff > 0);

        if (!newton_unconverged || iter == newton_max_iterations)
          {
            break;
          }

        // 5. Solve linearized system for the increment Delta u
        if (solutions.num_levels() > 1)
          {
            solutions.mg_transfer_down(
              solve_context->get_dof_manager(),
              solve_context->get_user_inputs().spatial_discretization.global_refinement,
              false);
          }

        newton_update = 0.0;
        total_lin_iters +=
          do_active_set_linear_solve(newton_residual, constrained_lhs_op, newton_update);
        zero_active_dofs(newton_update);
        newton_update.update_ghost_values();

        // 6. Backtracking line search keeping active DoFs pinned
        old_solution = solutions.get_solution_full_vector();
        number alpha = base_step_length;

        for (unsigned int ls = 0; ls < 5; ++ls)
          {
            solutions.get_solution_full_vector().zero_out_ghost_values();
            solutions.get_solution_full_vector() = old_solution;
            solutions.get_solution_full_vector().add(alpha, newton_update);
            solutions.apply_constraints();
            solutions.update_ghosts();

            // First iteration or active-set switch: accept step without line search
            if (iter == 0 || active_set_diff > 0)
              {
                break;
              }

            trial_residual.zero_out_ghost_values();
            rhs_op.compute_operator(trial_residual);
            trial_residual.update_ghost_values();
            zero_active_dofs(trial_residual);

            if (trial_residual.l2_norm() < l2_norm)
              {
                break;
              }

            alpha *= 0.5;
          }

        iter++;
      }

    if (solve_context->get_user_inputs().output_parameters.should_output(
          solve_context->get_simulation_timer().get_increment()))
      {
        const auto global_active = dealii::Utilities::MPI::sum(
          (active_set & locally_owned_dofs).n_elements(), MPI_COMM_WORLD);

        Logger::instance() << " Active-Set Newton final residual : "
                           << l2_norm / normalization_value()
                           << " | Active DoFs: " << global_active
                           << " | Newton steps: " << iter
                           << " | Total linear steps: " << total_lin_iters
                           << "\n"
                           << std::flush;
      }

    if (iter >= newton_max_iterations)
      {
        Logger::instance() << "[Increment "
                           << solve_context->get_simulation_timer().get_increment()
                           << "] "
                           << "Warning: Active-Set Newton solver did not converge before "
                           << newton_max_iterations << " iterations.\n\n";
      }
  }

protected:
  /**
   * @brief User-defined upper obstacle function psi(p) at support point p.
   * Override or bind to a dealii::Function<dim> for your specific problem.
   */
  virtual number
  evaluate_obstacle(const dealii::Point<dim> &p) const
  {
    (void)p;
    return static_cast<number>(0.0);
  }

private:
  unsigned int constrained_block;
  number       penalty_c;

  BlockVector<number> newton_update;
  BlockVector<number> old_solution;
  BlockVector<number> trial_residual;
  BlockVector<number> diag_mass_matrix;

  dealii::IndexSet locally_owned_dofs;
  dealii::IndexSet active_set;

  dealii::SolverControl                       active_set_solver_control;
  dealii::SolverSelector<BlockVector<number>> active_set_lin_solver;

  std::shared_ptr<dealii::DiagonalMatrix<BlockVector<number>>> active_set_inv_diag =
    std::make_shared<dealii::DiagonalMatrix<BlockVector<number>>>();

  [[nodiscard]] unsigned int
  constrained_field_index() const
  {
    return solutions.get_block_to_global_index()[constrained_block];
  }

  void
  initialize_active_set_linear_solver()
  {
    const auto &lp = lin_params();
    active_set_solver_control.set_max_steps(lp.max_iterations);
    active_set_solver_control.set_tolerance(lp.tolerance * normalization_value());

    const typename dealii::SolverRichardson<BlockVector<number>>::AdditionalData
      local_richardson_parameters(lp.richardson_parameters.omega,
                                  lp.richardson_parameters.use_preconditioned_residual);
    const typename dealii::SolverBicgstab<BlockVector<number>>::AdditionalData
      local_bicgstab_parameters(lp.bicgstab_parameters.exact_residual,
                                lp.bicgstab_parameters.breakdown);
    const typename dealii::SolverGMRES<BlockVector<number>>::AdditionalData
      local_gmres_parameters(lp.gmres_parameters.max_basis_size,
                             lp.gmres_parameters.right_preconditioning,
                             lp.gmres_parameters.use_default_residual,
                             lp.gmres_parameters.force_re_orthogonalization,
                             lp.gmres_parameters.batched_mode,
                             lp.gmres_parameters.orthogonalization_strategy);
    const typename dealii::SolverFGMRES<BlockVector<number>>::AdditionalData
      local_fgmres_parameters(lp.gmres_parameters.max_basis_size,
                              lp.gmres_parameters.orthogonalization_strategy);

    active_set_lin_solver.set_data(local_richardson_parameters);
    active_set_lin_solver.set_data(local_bicgstab_parameters);
    active_set_lin_solver.set_data(local_gmres_parameters);
    active_set_lin_solver.set_data(local_fgmres_parameters);

    active_set_lin_solver.select(lp.solver_type);
    active_set_lin_solver.set_control(active_set_solver_control);
  }

  int
  do_active_set_linear_solve(
    BlockVector<number>                           &b_vector,
    ActiveSetOperatorWrapper<dim, degree, number> &constrained_lhs_op,
    BlockVector<number>                           &x_vector)
  {
    // SolverSelector expects non-ghosted vectors
    b_vector.zero_out_ghost_values();
    x_vector.zero_out_ghost_values();

    try
      {
        if (lin_params().preconditioner == Chebyshev)
          {
            lhs_operator.reinit_matrix_diagonal();
            lhs_operator.eval_matrix_diagonal();

            // Copy inverse diagonal and set active-set rows to 1.0 to match
            // ActiveSetOperatorWrapper's identity rows.
            active_set_inv_diag->get_vector() =
              lhs_operator.get_matrix_diagonal_inverse()->get_vector();
            active_set_inv_diag->get_vector().zero_out_ghost_values();
            auto &diag_block =
              active_set_inv_diag->get_vector().block(constrained_block);
            for (const auto idx : active_set)
              {
                if (diag_block.in_local_range(idx))
                  {
                    diag_block[idx] = 1.0;
                  }
              }

            const auto &chebyshev_params = lin_params().chebyshev_parameters;
            typename ActiveSetPreconditionChebyshev::AdditionalData precond_data;
            precond_data.degree              = chebyshev_params.degree;
            precond_data.smoothing_range     = chebyshev_params.smoothing_range;
            precond_data.eig_cg_n_iterations = chebyshev_params.eig_cg_n_iterations;
            precond_data.preconditioner      = active_set_inv_diag;

            ActiveSetPreconditionChebyshev precond_chebyshev;
            precond_chebyshev.initialize(constrained_lhs_op, precond_data);

            active_set_lin_solver.solve(constrained_lhs_op,
                                        x_vector,
                                        b_vector,
                                        precond_chebyshev);
          }
        else
          {
            active_set_lin_solver.solve(constrained_lhs_op,
                                        x_vector,
                                        b_vector,
                                        dealii::PreconditionIdentity());
          }
      }
    catch (dealii::SolverControl::NoConvergence &)
      {
        Logger::instance()
          << "[Increment " << solve_context->get_simulation_timer().get_increment()
          << "] "
          << "Warning: active-set linear solver did not converge before "
          << lin_params().max_iterations << " iterations.\n";
      }

    return active_set_solver_control.last_step();
  }

  void
  setup_active_set_structures()
  {
    const auto &sol = solutions.get_solution_full_vector();
    newton_update.reinit(sol);
    old_solution.reinit(sol);
    trial_residual.reinit(sol);
    diag_mass_matrix.reinit(sol);

    locally_owned_dofs = sol.block(constrained_block).locally_owned_elements();
    active_set.clear();
    active_set.set_size(sol.block(constrained_block).size());

    compute_mass_matrix_diagonal();
  }

  /**
   * @brief Assemble diagonal lumped mass weights M_ii using QGaussLobatto
   * quadrature so support points coincide with quadrature points (step-42).
   */
  void
  compute_mass_matrix_diagonal()
  {
    const unsigned int field_index = constrained_field_index();
    const auto        &dof_handler =
      solve_context->get_dof_manager().get_field_dof_handler(field_index);
    const auto &fe = dof_handler.get_fe();

    // Use QGaussLobatto<dim> for a volume obstacle (e.g., phase-field bounds),
    // or QGaussLobatto<dim - 1> on boundary_id == 1 for surface contact (step-42).
    const dealii::QGaussLobatto<dim> quadrature(degree + 1);
    dealii::FEValues<dim>            fe_values(fe,
                                    quadrature,
                                    dealii::update_values |
                                      dealii::update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature.size();
    std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);

    diag_mass_matrix = 0.0;
    diag_mass_matrix.zero_out_ghost_values();

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (!cell->is_locally_owned())
          {
            continue;
          }

        fe_values.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                const number phi_i = fe_values.shape_value(i, q);
                diag_mass_matrix.block(constrained_block)[local_dof_indices[i]] +=
                  phi_i * phi_i * fe_values.JxW(q);
              }
          }
      }

    diag_mass_matrix.compress(dealii::VectorOperation::add);
    diag_mass_matrix.update_ghost_values();
  }

  /**
   * @brief Evaluate complementarity condition:
   *   lambda_i / M_ii + c * (u_i - psi_i) > 0
   * and project active DoFs onto u_i = psi_i.
   */
  void
  update_active_set_and_project(const BlockVector<number> &uncondensed_residual)
  {
    const unsigned int field_index = constrained_field_index();
    const auto        &dof_handler =
      solve_context->get_dof_manager().get_field_dof_handler(field_index);
    const auto &base_constraints =
      lhs_operator.get_matrix_free()->get_affine_constraints(field_index);
    const auto &fe = dof_handler.get_fe();

    const dealii::Quadrature<dim> support_quadrature(fe.get_unit_support_points());
    dealii::FEValues<dim>         fe_values(fe,
                                    support_quadrature,
                                    dealii::update_quadrature_points);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
    std::unordered_set<dealii::types::global_dof_index> visited_dofs;

    auto &sol_block = solutions.get_solution_full_vector().block(constrained_block);
    const auto &res_block  = uncondensed_residual.block(constrained_block);
    const auto &mass_block = diag_mass_matrix.block(constrained_block);

    active_set.clear();
    sol_block.zero_out_ghost_values();

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (!cell->is_locally_owned())
          {
            continue;
          }

        fe_values.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const auto dof_idx = local_dof_indices[i];

            if (!locally_owned_dofs.is_element(dof_idx) ||
                base_constraints.is_constrained(dof_idx) ||
                !visited_dofs.insert(dof_idx).second)
              {
                continue;
              }

            const number psi_i    = evaluate_obstacle(fe_values.quadrature_point(i));
            const number u_i      = sol_block[dof_idx];
            const number lambda_i = res_block[dof_idx];
            const number m_ii     = mass_block[dof_idx];

            if (m_ii > 0.0 && (lambda_i / m_ii + penalty_c * (u_i - psi_i) > 0.0))
              {
                active_set.add_index(dof_idx);
                sol_block[dof_idx] = psi_i;
              }
          }
      }

    sol_block.compress(dealii::VectorOperation::insert);
    solutions.update_ghosts();
  }

  void
  zero_active_dofs(BlockVector<number> &vec) const
  {
    vec.zero_out_ghost_values();
    auto &block = vec.block(constrained_block);
    for (const auto idx : active_set)
      {
        if (block.in_local_range(idx))
          {
            block[idx] = 0.0;
          }
      }
    vec.update_ghost_values();
  }

  [[nodiscard]] const NonlinearSolverParameters &
  newton_params() const
  {
    return solve_block.nonlinear_solver_parameters;
  }

  [[nodiscard]] const LinearSolverParameters &
  lin_params() const
  {
    return solve_block.linear_solver_parameters;
  }
};

PRISMS_PF_END_NAMESPACE