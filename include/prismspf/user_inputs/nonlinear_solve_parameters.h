// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/types.h>

#include <prismspf/user_inputs/linear_solve_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that stores relevant nonlinear solve information of a certain field
 */
struct NonlinearSolverParameters
{
  // Nonlinear step length
  double step_length = 1.0;

  // Max number of iterations for the nonlinear solve
  unsigned int max_iterations = Defaults::iterations;

  // Tolerance value for the nonlinear solve
  double tolerance_value = Defaults::tolerance;
};

/**
 * @brief Struct that holds nonlinear solver parameters.
 */
struct NonlinearSolveParameterSet
{
  /**
   * @brief Validate parameters.
   */
  void
  validate();

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler,
                     unsigned int              max_criteria = 5) const;

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  void
  assign_parameters(dealii::ParameterHandler &parameter_handler,
                    unsigned int              max_criteria = 5);

  // Map of nonlinear solve parameters for fields that require them
  std::map<Types::Index, NonlinearSolverParameters> newton_solvers;
};

inline void
NonlinearSolveParameterSet::validate()
{
  for (const auto &[index, nonlinear_solver_parameters] : newton_solvers)
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
NonlinearSolveParameterSet::declare_parameters(
  dealii::ParameterHandler &parameter_handler,
  unsigned int              max_criteria) const
{
  for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
    { // For nonlinear solves
      std::string subsection_text =
        "newton solver parameters: " + std::to_string(criterion_id);
      parameter_handler.enter_subsection(subsection_text);
      {
        parameter_handler.declare_entry(
          "solver_ids",
          "",
          dealii::Patterns::Anything(),
          "The ids of the solvers that will use these settings.");
        parameter_handler.declare_entry("max iterations",
                                        "100",
                                        dealii::Patterns::Integer(1, INT_MAX),
                                        "The maximum number of nonlinear solver "
                                        "iterations before the loop is stopped.");
        /*         parameter_handler.declare_entry(
                  "tolerance type",
                  "AbsoluteResidual",
                  dealii::Patterns::Selection("AbsoluteResidual|RMSEPerField|IntegratedPerField|"
                                              "RMSETotal|IntegratedTotal"),
                  "The tolerance type for the nonlinear solver."); */
        parameter_handler.declare_entry(
          "tolerance value",
          "1.0e-10",
          dealii::Patterns::Double(DBL_MIN, DBL_MAX),
          "The value of for the nonlinear solver tolerance.");
        /*         parameter_handler.declare_entry(
                  "use backtracking line search",
                  "true",
                  dealii::Patterns::Bool(),
                  "Whether to use a backtracking line-search to find the best "
                  "choice of the damping coefficient.");
                parameter_handler.declare_entry(
                  "step size modifier",
                  "0.5",
                  dealii::Patterns::Double(0.0, 1.0),
                  "The constant that determines how much the step size decreases "
                  "per backtrack. The 'tau' parameter.");
                parameter_handler.declare_entry(
                  "residual decrease coefficient",
                  "0.5",
                  dealii::Patterns::Double(0.0, 1.0),
                  "The constant that determines how much the residual must "
                  "decrease to be accepted as sufficient. The 'c' parameter."); */
        parameter_handler.declare_entry(
          "step size",
          "1.0",
          dealii::Patterns::Double(0.0),
          "The constant damping value to be used if the backtrace "
          "line-search approach isn't used.");
        parameter_handler.declare_alias("tolerance value", "tolerance");
        parameter_handler.declare_alias("solver_ids", "solve blocks");
        parameter_handler.declare_alias("solver_ids", "solve_blocks");
        parameter_handler.declare_alias("solver_ids", "solve block ids");
        parameter_handler.declare_alias("solver_ids", "solve_block_ids");
      }
      parameter_handler.leave_subsection();
    }
}

inline void
NonlinearSolveParameterSet::assign_parameters(dealii::ParameterHandler &parameter_handler,
                                              unsigned int              max_criteria)
{
  for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
    {
      // For newton solves
      std::string subsection_text =
        "newton solver parameters: " + std::to_string(criterion_id);
      parameter_handler.enter_subsection(subsection_text);
      {
        std::vector<int> solver_ids = dealii::Utilities::string_to_int(
          dealii::Utilities::split_string_list(parameter_handler.get("solver_ids")));
        NonlinearSolverParameters nonlinear_solver_parameters;
        nonlinear_solver_parameters.max_iterations =
          static_cast<unsigned int>(parameter_handler.get_integer("max iterations"));
        nonlinear_solver_parameters.step_length =
          parameter_handler.get_double("step size");
        nonlinear_solver_parameters.tolerance_value =
          parameter_handler.get_double("tolerance value");
        for (auto solver_id : solver_ids)
          {
            newton_solvers[static_cast<unsigned int>(solver_id)] =
              nonlinear_solver_parameters;
          }
      }
      parameter_handler.leave_subsection();
    }
}

PRISMS_PF_END_NAMESPACE
