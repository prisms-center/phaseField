// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_selector.h>

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
  // Solver type. richardson|cg|bicgstab|gmres|fgmres|minres
  std::string solver_type = "cg";

  // Solver tolerance
  double tolerance = Defaults::tolerance;

  // Solver tolerance type
  SolverToleranceType tolerance_type = SolverToleranceType::RMSEPerField;

  // Max number of iterations for the linear solve
  unsigned int max_iterations = Defaults::iterations;

  // Preconditioner
  PreconditionerType preconditioner = PreconditionerType::None;

  dealii::PreconditionChebyshev<>::AdditionalData chebyshev_parameters;

  dealii::SolverRichardson<>::AdditionalData richardson_parameters;
  dealii::SolverBicgstab<>::AdditionalData   bicgstab_parameters;
  dealii::SolverGMRES<>::AdditionalData      gmres_parameters;

  // MinRes and CG do not have additional parameters, so we do not need to store them here
  // dealii::SolverMinRes<>::AdditionalData minres_parameters;
  // dealii::SolverCG<>::AdditionalData     cg_parameters;
  // FGMRES additional parameters are a subset of GMRES parameters, so we can just use
  // gmres_parameters for both solvers
  // dealii::SolverFGMRES<>::AdditionalData     fgmres_parameters;

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
              // ConditionalOStreams::pout_summary();
            }
        }

      ConditionalOStreams::pout_summary() << "\n" << std::flush;
    }
}

PRISMS_PF_END_NAMESPACE
