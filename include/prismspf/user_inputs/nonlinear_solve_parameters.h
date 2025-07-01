// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/types.h>

#include <prismspf/user_inputs/linear_solve_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Struct that stores relevant nonlinear solve information of a certain field
 */
struct NonlinearSolverParameters
{
public:
  // Nonlinear step length
  mutable double step_length = 1.0;

  // Max number of iterations for the nonlinear solve
  unsigned int max_iterations = Defaults::iterations;

  // Tolerance value for the nonlinear solve
  double tolerance_value = Defaults::tolerance;
};

/**
 * \brief Struct that holds nonlinear solver parameters.
 */
struct NonlinearSolveParameterSet
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
   * \brief Clear the nonlinear solve parameters.
   */
  void
  clear()
  {
    nonlinear_solve.clear();
  }

  /**
   * \brief Set the nonlinear solve parameters for a field index.
   */
  void
  set_nonlinear_solve_parameters(Types::Index                     index,
                                 const NonlinearSolverParameters &parameters)
  {
    nonlinear_solve[index] = parameters;
  }

  /**
   * \brief Get the nonlinear solve parameters for a field index.
   */
  [[nodiscard]] const NonlinearSolverParameters &
  get_nonlinear_solve_parameters(Types::Index index) const
  {
    return nonlinear_solve.at(index);
  }

private:
  // Map of nonlinear solve parameters for fields that require them
  std::map<Types::Index, NonlinearSolverParameters> nonlinear_solve;
};

inline void
NonlinearSolveParameterSet::postprocess_and_validate()
{
  for (const auto &[index, nonlinear_solver_parameters] : nonlinear_solve)
    {
      AssertThrow(
        nonlinear_solver_parameters.step_length > 0.0 &&
          nonlinear_solver_parameters.step_length <= 1.0,
        dealii::ExcMessage(
          "Step length must be greater than 0.0 and less than or equal to 1.0"));

      AssertThrow(nonlinear_solver_parameters.tolerance_value > 0,
                  dealii::ExcMessage("Tolerance must be greater than 0.0"));
    }
}

inline void
NonlinearSolveParameterSet::print_parameter_summary() const
{
  if (!nonlinear_solve.empty())
    {
      ConditionalOStreams::pout_summary()
        << "================================================\n"
        << "  Nonlinear Solve Parameters\n"
        << "================================================\n";

      for (const auto &[index, nonlinear_solver_parameters] : nonlinear_solve)
        {
          ConditionalOStreams::pout_summary()
            << "Index: " << index << "\n"
            << "  Max iterations: " << nonlinear_solver_parameters.max_iterations << "\n"
            << "  Step length: " << nonlinear_solver_parameters.step_length << "\n";
        }

      ConditionalOStreams::pout_summary() << "\n" << std::flush;
    }
}

PRISMS_PF_END_NAMESPACE