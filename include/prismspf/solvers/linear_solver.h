// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_selector.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/invm_manager.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/mf_operator.h>
#include <prismspf/solvers/solver_base.h>

#include <prismspf/user_inputs/linear_solve_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolveContext;

/**
 * @brief This class handles the explicit solves of all explicit fields
 * @todo Separate implementation from header.
 */
template <unsigned int dim, unsigned int degree, typename number>
class LinearSolver : public SolverBase<dim, degree, number>
{
protected:
  using SolverBase<dim, degree, number>::solutions;
  using SolverBase<dim, degree, number>::solve_context;
  using SolverBase<dim, degree, number>::solve_block;
  using PreconditionChebyshev =
    dealii::PreconditionChebyshev<MFOperator<dim, degree, number>,
                                  BlockVector<number>,
                                  dealii::DiagonalMatrix<BlockVector<number>>>;

public:
  /**
   * @brief Constructor.
   * @pre Solve context has initialized members.
   */
  LinearSolver(SolveBlock                               _solve_block,
               const SolveContext<dim, degree, number> &_solve_context)
    : SolverBase<dim, degree, number>(_solve_block, _solve_context)
    , rhs_operator(solve_context->get_pde_operator(),
                   &PDEOperatorBase<dim, degree, number>::compute_rhs,
                   solve_context->get_field_attributes(),
                   solve_context->get_solution_indexer(),
                   solve_context->get_matrix_free_manager(),
                   solve_block.dependencies_rhs,
                   solve_context->get_simulation_timer())
    , lhs_operator(solve_context->get_pde_operator(),
                   &PDEOperatorBase<dim, degree, number>::compute_lhs,
                   solve_context->get_field_attributes(),
                   solve_context->get_solution_indexer(),
                   solve_context->get_matrix_free_manager(),
                   solve_block.dependencies_lhs,
                   solve_context->get_simulation_timer())
  {}

  /**
   * @brief Initialize the solver.
   */
  void
  init(const std::list<DependencyMap> &all_dependeny_sets) override
  {
    SolverBase<dim, degree, number>::init(all_dependeny_sets);
    rhs_vector.reinit(solutions.get_solution_full_vector(0));

    // Initialize rhs_operator
    rhs_operator.initialize(solutions);
    rhs_operator.set_scaling_diagonal(lin_params().tolerance_type != AbsoluteResidual,
                                      solve_context->get_invm_manager().get_invm_sqrt(
                                        solve_context->get_field_attributes(),
                                        solve_block.field_indices,
                                        0));
    // Initialize lhs_operator
    lhs_operator.initialize(solutions);
    lhs_operator.set_scaling_diagonal(lin_params().tolerance_type != AbsoluteResidual,
                                      solve_context->get_invm_manager().get_invm_sqrt(
                                        solve_context->get_field_attributes(),
                                        solve_block.field_indices,
                                        0));

    linear_solver_control.set_max_steps(lin_params().max_iterations);
    linear_solver_control.set_tolerance(lin_params().tolerance * normalization_value());
    initialize_solver();
    inhomogeneous_values.reinit(solutions.get_solution_full_vector(0));
    solutions.apply_constraints(inhomogeneous_values, 0);
    inhomogeneous_rhs.reinit(solutions.get_solution_full_vector(0));
    initialize_preconditioner();
  }

  /**
   * @brief Reinitialize the solver.
   */
  void
  reinit() override
  {
    SolverBase<dim, degree, number>::reinit();
    rhs_vector.reinit(solutions.get_solution_full_vector(0));

    inhomogeneous_values.reinit(solutions.get_solution_full_vector(0));
    solutions.apply_constraints(inhomogeneous_values, 0);
    inhomogeneous_rhs.reinit(solutions.get_solution_full_vector(0));
  }

  /**
   * @brief Solve for a single update step.
   */
  void
  solve_impl() override
  {
    // Zero out the ghosts
    Timer::start_section("Zero ghosts");
    solutions.zero_out_ghosts(0);
    Timer::end_section("Zero ghosts");

    // Set up rhs vector
    rhs_operator.compute_operator(rhs_vector);

    // Note 1. Use the previous result of the linear solve without nonzero dirichlet
    // as the initial guess in the next increment. See Note 2. `inhomogeneous_rhs` is
    // not actually what it is being used as here, we just don't want to allocate a
    // whole new vector for this purpose
    solutions.get_solution_full_vector(0).swap(inhomogeneous_rhs);
    // Get the homogeneous rhs
    lhs_operator.read_plain = true;
    lhs_operator.compute_operator(inhomogeneous_rhs, inhomogeneous_values);
    lhs_operator.read_plain = false;
    rhs_vector -= inhomogeneous_rhs;

    // Linear solve
    do_linear_solve(rhs_vector, lhs_operator, solutions.get_solution_full_vector(0));

    // Note 2. Make a copy of the solution to use as the initial guess in the next
    // increment. See Note 1. `inhomogeneous_rhs` is not actually what it is being
    // used as here, we just don't want to allocate a whole new vector for this
    // purpose
    inhomogeneous_rhs = solutions.get_solution_full_vector(0);
    // Add back in nonzero dirichlet conditions
    solutions.get_solution_full_vector(0) += inhomogeneous_values;

    // Apply constraints
    solutions.apply_constraints(0);

    // Update the ghosts
    Timer::start_section("Update ghosts");
    solutions.update_ghosts(0);
    Timer::end_section("Update ghosts");
  }

  int
  do_linear_solve(BlockVector<number>             &b_vector,
                  MFOperator<dim, degree, number> &lhs_matrix,
                  BlockVector<number>             &x_vector)
  {
    // Linear solve
    try
      {
        if (lin_params().preconditioner == None)
          {
            lin_solver.solve(lhs_matrix,
                             x_vector,
                             b_vector,
                             dealii::PreconditionIdentity());
          }
        if (lin_params().preconditioner == Chebyshev)
          {
            lhs_matrix.reinit_matrix_diagonal(x_vector);
            lhs_matrix.eval_matrix_diagonal();

            lin_solver.solve(lhs_matrix, x_vector, b_vector, precond_chebyshev);
          }
      }
    catch (dealii::SolverControl::NoConvergence &exc)
      {
        ConditionalOStreams::pout_base()
          << "[Increment " << solve_context->get_simulation_timer().get_increment()
          << "] "
          << "Warning: linear solver did not converge as per set tolerances before "
          << lin_params().max_iterations << " iterations.\n";
      }
    if (solve_context->get_user_inputs().output_parameters.should_output(
          solve_context->get_simulation_timer().get_increment()))
      {
        ConditionalOStreams::pout_summary()
          << " Linear solve final residual : "
          << linear_solver_control.last_value() / normalization_value()
          << " Linear steps: " << linear_solver_control.last_step() << "\n"
          << std::flush;
      }
    return linear_solver_control.last_step();
  }

protected:
  /**
   * @brief Matrix free operators
   */
  MFOperator<dim, degree, number> rhs_operator;

  MFOperator<dim, degree, number> lhs_operator;
  BlockVector<number>             rhs_vector;

  double
  normalization_value()
  {
    SolverToleranceType type  = lin_params().tolerance_type;
    double              value = 1.0;
    if (type == RMSEPerField || type == RMSETotal)
      {
        value *= std::sqrt(solve_context->get_triangulation_manager().get_volume());
      }
    if (type == RMSEPerField || type == IntegratedPerField)
      {
        value *= std::sqrt(double(solve_block.field_indices.size()));
      }
    return value;
  }

  void
  initialize_preconditioner()
  {
    if (lin_params().preconditioner == None)
      {
        void(0); // do nothing
      }
    if (lin_params().preconditioner == Chebyshev)
      {
        initialize_chebyshev();
      }
  }

  void
  initialize_chebyshev()
  {
    const auto &chebyshev_params = lin_params().chebyshev_parameters;
    precond_data.degree          = chebyshev_params.degree;
    precond_data.smoothing_range = chebyshev_params.smoothing_range; // ≈ λ_min / λ_max
    precond_data.eig_cg_n_iterations = chebyshev_params.eig_cg_n_iterations;

    lhs_operator.reinit_matrix_diagonal(solutions.get_solution_full_vector(0));
    precond_data.preconditioner = lhs_operator.get_matrix_diagonal_inverse();

    precond_chebyshev.initialize(lhs_operator, precond_data);
  }

  void
  initialize_solver()
  {
    const auto &richardson_parameters = lin_params().richardson_parameters;
    const auto &bicgstab_parameters   = lin_params().bicgstab_parameters;
    const auto &gmres_parameters      = lin_params().gmres_parameters;
    const typename dealii::SolverRichardson<BlockVector<number>>::AdditionalData
      local_richardson_parameters(richardson_parameters.omega,
                                  richardson_parameters.use_preconditioned_residual);
    const typename dealii::SolverBicgstab<BlockVector<number>>::AdditionalData
      local_bicgstab_parameters(bicgstab_parameters.exact_residual,
                                bicgstab_parameters.breakdown);
    const typename dealii::SolverGMRES<BlockVector<number>>::AdditionalData
      local_gmres_parameters(gmres_parameters.max_basis_size,
                             gmres_parameters.right_preconditioning,
                             gmres_parameters.use_default_residual,
                             gmres_parameters.force_re_orthogonalization,
                             gmres_parameters.batched_mode,
                             gmres_parameters.orthogonalization_strategy);
    const typename dealii::SolverFGMRES<BlockVector<number>>::AdditionalData
      local_fgmres_parameters(gmres_parameters.max_basis_size,
                              gmres_parameters.orthogonalization_strategy);

    lin_solver.set_data(local_richardson_parameters);
    lin_solver.set_data(local_bicgstab_parameters);
    lin_solver.set_data(local_gmres_parameters);
    lin_solver.set_data(local_fgmres_parameters);

    lin_solver.select(lin_params().solver_type);
    lin_solver.set_control(linear_solver_control);
  }

private:
  /**
   * @brief Linear solver parameters
   */
  const LinearSolverParameters &
  lin_params() const
  {
    return solve_block.linear_solver_parameters;
  }

  /**
   * @brief Solver control. Contains max iterations and tolerance.
   */
  dealii::SolverControl linear_solver_control;

  /**
   * @brief Solver. Can switch between different linear solvers.
   */
  dealii::SolverSelector<BlockVector<number>> lin_solver;

  /**
   * @brief Vector containing only the inhomogeneous constraints (namely, non-zero
   * Dirichlet values)
   */
  BlockVector<number> inhomogeneous_values;

  /**
   * @brief Result of the linear operator applied to the inhomogeneous values.
   */
  BlockVector<number> inhomogeneous_rhs;

  PreconditionChebyshev                 precond_chebyshev;
  PreconditionChebyshev::AdditionalData precond_data;
};

PRISMS_PF_END_NAMESPACE
