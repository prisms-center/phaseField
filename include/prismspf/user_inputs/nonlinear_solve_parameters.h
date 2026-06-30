// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/types.h>

#include <prismspf/user_inputs/linear_solve_parameters.h>
#include <prismspf/user_inputs/parameter_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that stores relevant nonlinear solve information of a certain field
 */
struct NonlinearSolverParameters : public ParameterBase
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int max_criteria = Numbers::max_subsections) const override
  {
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

    parameter_handler.declare_entry("tolerance value",
                                    "1.0e-10",
                                    dealii::Patterns::Double(0.0, DBL_MAX),
                                    "The value of for the nonlinear solver tolerance.");
    parameter_handler.declare_alias("tolerance value", "tolerance");

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
  };

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              max_criteria = Numbers::max_subsections) override
  {
    max_iterations = (unsigned int) (parameter_handler.get_integer("max iterations"));

    step_length = parameter_handler.get_double("step size");

    tolerance_value = parameter_handler.get_double("tolerance value");
  };

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override {
    // TODO: Add validation here
  };

  // Nonlinear step length
  double step_length = 1.0;

  // Max number of iterations for the nonlinear solve
  unsigned int max_iterations = 100;

  // Tolerance value for the nonlinear solve
  double tolerance_value = 1.0e-10;
};

/**
 * @brief Struct that holds nonlinear solver parameters.
 */
struct NonlinearSolveParameters : public ParameterBase
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int max_criteria = Numbers::max_subsections) const override
  {
    for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
      {
        std::string subsection_text =
          "newton solver parameters: " + std::to_string(criterion_id);
        parameter_handler.enter_subsection(subsection_text);
        {
          parameter_handler.declare_entry(
            "solver_ids",
            "",
            dealii::Patterns::List(dealii::Patterns::Anything(), 0, INT_MAX, ","),
            "The ids of the solvers that will use these settings.");
          declare_aliases(parameter_handler,
                          "solver_ids",
                          {"solve blocks",
                           "solve_blocks",
                           "solve block ids",
                           "solve_block_ids",
                           "solver ids"});

          NonlinearSolverParameters nonlinear_solver;
          nonlinear_solver.declare(parameter_handler);
        }
        parameter_handler.leave_subsection();
      }
  };

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              max_criteria = Numbers::max_subsections) override
  {
    for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
      {
        std::string subsection_text =
          "newton solver parameters: " + std::to_string(criterion_id);
        parameter_handler.enter_subsection(subsection_text);
        {
          std::vector<int> solver_ids = dealii::Utilities::string_to_int(
            dealii::Utilities::split_string_list(parameter_handler.get("solver_ids")));

          NonlinearSolverParameters nonlinear_solver;
          nonlinear_solver.assign(parameter_handler);

          for (auto solver_id : solver_ids)
            {
              nonlinear_solvers[(unsigned int) (solver_id)] = nonlinear_solver;
            }
        }
        parameter_handler.leave_subsection();
      }
  };

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override
  {
    for (const auto &[solver_id, nonlinear_solver] : nonlinear_solvers)
      {
        nonlinear_solver.validate(field_attributes, solve_blocks);
      }
  };

  // Map of nonlinear solve parameters for fields that require them
  std::map<Types::Index, NonlinearSolverParameters> nonlinear_solvers;
};

PRISMS_PF_END_NAMESPACE
