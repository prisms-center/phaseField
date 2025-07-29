// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/sequential_solver.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class handles the co-nonlinear solves of a several nonexplicit fields
 */
template <unsigned int dim, unsigned int degree, typename number>
class SequentialCoNonlinearSolver : public SequentialSolver<dim, degree, number>
{
public:
  /**
   * @brief Constructor.
   */
  explicit SequentialCoNonlinearSolver(const SolverContext<dim, degree> &_solver_context,
                                       Types::Index _solve_priority = 0)
    : SequentialSolver<dim, degree, number>(_solver_context,
                                            FieldSolveType::NonexplicitCononlinear,
                                            _solve_priority) {};

  /**
   * @brief Destructor.
   */
  ~SequentialCoNonlinearSolver() override = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialCoNonlinearSolver(const SequentialCoNonlinearSolver &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialCoNonlinearSolver &
  operator=(const SequentialCoNonlinearSolver &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialCoNonlinearSolver(SequentialCoNonlinearSolver &&solver_base) noexcept =
    delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialCoNonlinearSolver &
  operator=(SequentialCoNonlinearSolver &&solver_base) noexcept = delete;

  /**
   * @brief Initialize the solver.
   */
  void
  init() override
  {
    // Call the base class init
    this->SequentialSolver<dim, degree, number>::init();

    // If the solver is empty we can just return early.
    if (this->solver_is_empty())
      {
        return;
      }

    // Init the linear and auxiliary solvers
    for (const auto &[index, variable] : this->get_subset_attributes())
      {
        if (variable.get_pde_type() == PDEType::Auxiliary)
          {
            this->init_explicit_solver(variable);
          }
        else if (variable.get_pde_type() == PDEType::ImplicitTimeDependent ||
                 variable.get_pde_type() == PDEType::TimeIndependent)
          {
            this->init_linear_solver(variable);
          }
        else
          {
            AssertThrow(false, UnreachableCode());
          }
      }
  };

  /**
   * @brief Reinitialize the solver.
   */
  void
  reinit() override
  {
    // Call the base class reinit
    this->SequentialSolver<dim, degree, number>::reinit();

    // If the solver is empty we can just return early.
    if (this->solver_is_empty())
      {
        return;
      }

    // Reinit the linear and auxiliary solvers
    for (const auto &[index, variable] : this->get_subset_attributes())
      {
        if (variable.get_pde_type() == PDEType::Auxiliary)
          {
            this->reinit_explicit_solver(variable);
          }
        else if (variable.get_pde_type() == PDEType::ImplicitTimeDependent ||
                 variable.get_pde_type() == PDEType::TimeIndependent)
          {
            this->reinit_linear_solver(variable);
          }
        else
          {
            AssertThrow(false, UnreachableCode());
          }
      }
  };

  /**
   * @brief Solve for a single update step.
   */
  void
  solve() override
  {
    // Call the base class solve
    this->SequentialSolver<dim, degree, number>::solve();

    // If the solver is empty we can just return early.
    if (this->solver_is_empty())
      {
        return;
      }

    // Set the convergence bool and iteration counter
    bool         unconverged = true;
    unsigned int iteration   = 0;

    while (unconverged)
      {
        if (this->get_user_inputs().get_output_parameters().should_output(
              this->get_user_inputs().get_temporal_discretization().get_increment()))
          {
            ConditionalOStreams::pout_summary()
              << "Nonlinear solver step: " << iteration << "\n";
          }

        // Assume the solve is converged, unless proven otherwise
        unconverged = false;

        for (const auto &[index, variable] : this->get_subset_attributes())
          {
            // Set the step length
            const double step_length = this->get_user_inputs()
                                         .get_nonlinear_solve_parameters()
                                         .get_nonlinear_solve_parameters(index)
                                         .step_length;

            // Set the norm of the newton update
            double newton_update_norm = 0.0;

            if (variable.get_pde_type() == PDEType::Auxiliary)
              {
                this->solve_explicit_solver(variable);
              }
            else if (variable.get_pde_type() == PDEType::ImplicitTimeDependent ||
                     variable.get_pde_type() == PDEType::TimeIndependent)
              {
                // Perform the linear solve with the step length
                newton_update_norm = this->solve_linear_solver(variable, step_length);
              }
            else
              {
                AssertThrow(false, UnreachableCode());
              }

            // Check the convergence of the nonlinear solve
            if (this->get_user_inputs().get_output_parameters().should_output(
                  this->get_user_inputs().get_temporal_discretization().get_increment()))
              {
                ConditionalOStreams::pout_summary()
                  << "  field: " << index << " Newton update norm: " << newton_update_norm
                  << "\n"
                  << std::flush;
              }

            if (newton_update_norm > this->get_user_inputs()
                                       .get_nonlinear_solve_parameters()
                                       .get_nonlinear_solve_parameters(index)
                                       .tolerance_value)
              {
                unconverged = true;
              }

            // Check if the maximum number of iterations has been reached
            if (iteration >= this->get_user_inputs()
                               .get_nonlinear_solve_parameters()
                               .get_nonlinear_solve_parameters(index)
                               .max_iterations)
              {
                unconverged = false;
                ConditionalOStreams::pout_base() << "Warning: nonlinear solver did not "
                                                    "converge as per set tolerances.\n\n"
                                                 << std::flush;
              }
          }

        // Update the iteration counter
        iteration++;
      }

    // Update the solutions
    for (const auto &[index, variable] : this->get_subset_attributes())
      {
        this->get_solution_handler().update(this->get_field_solve_type(),
                                            this->get_solve_block(),
                                            index);
      }

    // Update the ghosts
    Timer::start_section("Update ghosts");
    this->get_solution_handler().update_ghosts();
    Timer::end_section("Update ghosts");
  };

  /**
   * @brief Print information about the solver to summary.log.
   */
  void
  print() override
  {
    // Print the base class information
    this->SequentialSolver<dim, degree, number>::print();
  }
};

PRISMS_PF_END_NAMESPACE
