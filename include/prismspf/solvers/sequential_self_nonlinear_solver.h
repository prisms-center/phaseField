// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/sequential_solver.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class handles the self-nonlinear solves of a single nonexplicit field
 */
template <unsigned int dim, unsigned int degree, typename number>
class SequentialSelfNonlinearSolver : public SequentialSolver<dim, degree, number>
{
public:
  /**
   * @brief Constructor.
   */
  explicit SequentialSelfNonlinearSolver(
    const SolverContext<dim, degree> &_solver_context,
    Types::Index                      _solve_priority = 0)
    : SequentialSolver<dim, degree, number>(_solver_context,
                                            FieldSolveType::NonexplicitSelfnonlinear,
                                            _solve_priority) {};

  /**
   * @brief Destructor.
   */
  ~SequentialSelfNonlinearSolver() override = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialSelfNonlinearSolver(const SequentialSelfNonlinearSolver &solver_base) =
    delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialSelfNonlinearSolver &
  operator=(const SequentialSelfNonlinearSolver &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialSelfNonlinearSolver(SequentialSelfNonlinearSolver &&solver_base) noexcept =
    delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialSelfNonlinearSolver &
  operator=(SequentialSelfNonlinearSolver &&solver_base) noexcept = delete;

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

    // Init the linear solvers
    for (const auto &[index, variable] : this->get_subset_attributes())
      {
        this->init_linear_solver(variable);
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

    // Solve each field
    for (const auto &[index, variable] : this->get_subset_attributes())
      {
        // Set the convergence bool, iteration counter, and step length
        bool         unconverged = true;
        unsigned int iteration   = 0;
        const double step_length = this->get_user_inputs()
                                     .get_nonlinear_solve_parameters()
                                     .get_nonlinear_solve_parameters(index)
                                     .step_length;

        while (unconverged)
          {
            if (this->get_user_inputs().get_output_parameters().should_output(
                  this->get_user_inputs()
                    .get_temporal_discretization()
                    .get_current_increment()))
              {
                ConditionalOStreams::pout_summary()
                  << "Nonlinear solver step: " << iteration << "\n";
              }

            // Assume the solve is converged, unless proven otherwise
            unconverged = false;

            // Set the norm of the newton update
            double newton_update_norm = 0.0;

            // Perform the linear solve with the step length
            newton_update_norm = this->solve_linear_solver(variable, step_length);

            // Check the convergence of the nonlinear solve
            if (this->get_user_inputs().get_output_parameters().should_output(
                  this->get_user_inputs()
                    .get_temporal_discretization()
                    .get_current_increment()))
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

            // Update the iteration counter
            iteration++;
          }

        // Update the solutions
        this->get_solution_handler().update(this->get_field_solve_type(), index);

        // Update the ghosts
        Timer::start_section("Update ghosts");
        this->get_solution_handler().update_ghosts();
        Timer::end_section("Update ghosts");
      }
  };

  /**
   * @brief Print information about the solver to summary.log.
   */
  void
  print()
  {
    // Print the base class information
    this->SequentialSolver<dim, degree, number>::print();
  }
};

PRISMS_PF_END_NAMESPACE
