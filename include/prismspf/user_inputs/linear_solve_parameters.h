// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <map>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that stores relevant linear solve information of a certain field
 */
struct LinearSolverParameters
{
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
 * @brief Struct that holds linear solver parameters.
 */
struct LinearSolveParameters
{
  /**
   * @brief Postprocess and validate parameters.
   */
  void
  validate();

  /**
   * @brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  // Map of linear solve parameters for fields that require them
  std::map<unsigned int, LinearSolverParameters> linear_solvers;
};

inline void
LinearSolveParameters::validate()
{
  // Nothing to do here for now
}

inline void
LinearSolveParameters::print_parameter_summary() const
{
  if (!linear_solvers.empty())
    {
      ConditionalOStreams::pout_summary()
        << "================================================\n"
        << "  Linear Solve Parameters\n"
        << "================================================\n";

      for (const auto &[index, linear_solver_parameters] : linear_solvers)
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
