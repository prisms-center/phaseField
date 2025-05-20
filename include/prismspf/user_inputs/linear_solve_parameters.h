// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <map>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Struct that stores relevant linear solve information of a certain field
 */
struct LinearSolverParameters
{
public:
  // Solver tolerance
  double tolerance = Defaults::tolerance;

  // Solver tolerance type
  SolverToleranceType tolerance_type = SolverToleranceType::RelativeResidualChange;

  // Max number of iterations for the linear solve
  unsigned int max_iterations = Defaults::iterations;

  // Preconditioner
  PreconditionerType preconditioner = PreconditionerType::GMG;

  // Smoothing range for eigenvalues. This denotes the lower bound of eigenvalues that are
  // smoothed [1.2 λ^max / smoothing_range, 1.2 λ^max], where λ^max is the estimated
  // maximum eigenvalue. A choice between 5 and 20 is usually useful when the
  // preconditioner is used as a smoother in multigrid.
  double smoothing_range = Defaults::smoothing_range;

  // Polynomial degree for the Chebyshev smoother
  unsigned int smoother_degree = Defaults::smoother_degree;

  // Maximum number of CG iterations used to find the maximum eigenvalue
  unsigned int eig_cg_n_iterations = Defaults::eig_cg_n_iterations;

  // The minimum multigrid level
  unsigned int min_mg_level = 0;
};

/**
 * \brief Struct that holds linear solver parameters.
 */
struct LinearSolveParameters
{
public:
  /**
   * \brief Postprocess and validate parameters.
   */
  void
  postprocess_and_validate();

  /**
   * \brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  /**
   * \brief Whether we have any linear solve parameters.
   */
  [[nodiscard]] bool
  has_linear_solve_parameters() const
  {
    return !linear_solve.empty();
  }

  /**
   * \brief Whether we have linear solve parameters for a given field.
   */
  [[nodiscard]] bool
  has_linear_solve_parameters(unsigned int field_index) const
  {
    return linear_solve.contains(field_index);
  }

  /**
   * \brief Set the linear solve parameters for a given field.
   */
  void
  set_linear_solve_parameters(unsigned int                  field_index,
                              const LinearSolverParameters &linear_solver_parameters)
  {
    linear_solve[field_index] = linear_solver_parameters;
  }

  /**
   * \brief Return the linear solve parameters for a given field.
   */
  [[nodiscard]] const LinearSolverParameters &
  get_linear_solve_parameters(unsigned int field_index) const
  {
    AssertThrow(has_linear_solve_parameters(field_index),
                dealii::ExcMessage("No linear solve parameters found for field index " +
                                   std::to_string(field_index)));
    return linear_solve.at(field_index);
  }

  /**
   * \brief Return the linear solve parameters for a given field.
   */
  [[nodiscard]] const std::map<unsigned int, LinearSolverParameters> &
  get_linear_solve_parameters() const
  {
    return linear_solve;
  }

private:
  // Map of linear solve parameters for fields that require them
  std::map<unsigned int, LinearSolverParameters> linear_solve;
};

inline void
LinearSolveParameters::postprocess_and_validate()
{
  // Nothing to do here for now
}

inline void
LinearSolveParameters::print_parameter_summary() const
{
  if (!linear_solve.empty())
    {
      ConditionalOStreams::pout_summary()
        << "================================================\n"
        << "  Linear Solve Parameters\n"
        << "================================================\n";

      for (const auto &[index, linear_solver_parameters] : linear_solve)
        {
          ConditionalOStreams::pout_summary()
            << "Index: " << index << "\n"
            << "  Tolerance: " << linear_solver_parameters.tolerance << "\n"
            << "  Type: " << to_string(linear_solver_parameters.tolerance_type) << "\n"
            << "  Max iterations: " << linear_solver_parameters.max_iterations << "\n"
            << "  Preconditioner: " << to_string(linear_solver_parameters.preconditioner)
            << "\n";

          if (linear_solver_parameters.preconditioner == PreconditionerType::GMG)
            {
              ConditionalOStreams::pout_summary()
                << "  Smoothing range: " << linear_solver_parameters.smoothing_range
                << "\n"
                << "  Smoother degree: " << linear_solver_parameters.smoother_degree
                << "\n"
                << "  Max eigenvalue CG iterations: "
                << linear_solver_parameters.eig_cg_n_iterations << "\n"
                << "Min multigrid level: " << linear_solver_parameters.min_mg_level
                << "\n";
            }
        }

      ConditionalOStreams::pout_summary() << "\n" << std::flush;
    }
}

PRISMS_PF_END_NAMESPACE
