// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_block.h>

#include <prismspf/user_inputs/constraint_parameters.h>
#include <prismspf/user_inputs/io_parameters.h>
#include <prismspf/user_inputs/miscellaneous_parameters.h>
#include <prismspf/user_inputs/nucleation_parameters.h>
#include <prismspf/user_inputs/solve_parameters.h>
#include <prismspf/user_inputs/spatial_discretization.h>
#include <prismspf/user_inputs/temporal_discretization.h>
#include <prismspf/user_inputs/user_constants.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
struct UserInputParameters : public ParameterBase
{
  /**
   * @brief Constructor.
   */
  UserInputParameters() = default;

  /**
   * @brief Constructor.
   *
   * Read in user input parameters from file and load them into member variables.
   */
  explicit UserInputParameters(const std::string &file_name,
                               unsigned int       max_criteria = Numbers::max_subsections)
  {
    // user_constants is a special little princess that needs to know the file contents
    // before we declare the parameters.
    user_constants.file_name = file_name;
    dealii::ParameterHandler parameter_handler;
    predeclare(parameter_handler);
    parameter_handler.parse_input(file_name);
    preassign(parameter_handler);
    parameter_handler.clear();
    declare(parameter_handler, max_criteria);
    parameter_handler.parse_input(file_name);
    assign(parameter_handler, max_criteria);
  };

  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  predeclare(dealii::ParameterHandler &parameter_handler) const override
  {
    spatial_discretization.predeclare(parameter_handler);
  };

  /**
   * @brief Assign the parameters from file.
   */
  void
  preassign(dealii::ParameterHandler &parameter_handler) override
  {
    spatial_discretization.preassign(parameter_handler);
  };

  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int max_criteria = Numbers::max_subsections) const override
  {
    spatial_discretization.declare(parameter_handler, max_criteria);
    temporal_discretization.declare(parameter_handler, max_criteria);
    boundary_parameters.declare(parameter_handler, max_criteria);

    linear_solve_parameters.declare(parameter_handler, max_criteria);
    nonlinear_solve_parameters.declare(parameter_handler, max_criteria);

    output_parameters.declare(parameter_handler, max_criteria);
    restart_parameters.declare(parameter_handler, max_criteria);
    input_parameters.declare(parameter_handler, max_criteria);

    misc_parameters.declare(parameter_handler, max_criteria);
    nucleation_parameters.declare(parameter_handler, max_criteria);

    user_constants.declare_parameters(parameter_handler);
  };

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              max_criteria = Numbers::max_subsections) override
  {
    spatial_discretization.assign(parameter_handler, max_criteria);
    temporal_discretization.assign(parameter_handler, max_criteria);
    boundary_parameters.assign(parameter_handler, max_criteria);

    linear_solve_parameters.assign(parameter_handler, max_criteria);
    nonlinear_solve_parameters.assign(parameter_handler, max_criteria);

    const auto n_increments = temporal_discretization.n_increments;

    output_parameters.assign(parameter_handler, n_increments, max_criteria);
    restart_parameters.assign(parameter_handler, n_increments, max_criteria);
    input_parameters.assign(parameter_handler, max_criteria);

    misc_parameters.assign(parameter_handler, max_criteria);
    nucleation_parameters.assign(parameter_handler, max_criteria);

    user_constants.assign_parameters(parameter_handler);
  };

  /**
   * @brief Ensure that the parameters are compatible with a set of fields and solvers.
   *
   * This is an optional step.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override
  {
    spatial_discretization.validate(field_attributes, solve_blocks);
    temporal_discretization.validate(field_attributes, solve_blocks);
    boundary_parameters.validate(field_attributes, solve_blocks);

    linear_solve_parameters.validate(field_attributes, solve_blocks);
    nonlinear_solve_parameters.validate(field_attributes, solve_blocks);

    output_parameters.validate(field_attributes, solve_blocks);
    restart_parameters.validate(field_attributes, solve_blocks);
    input_parameters.validate(field_attributes, solve_blocks);

    misc_parameters.validate(field_attributes, solve_blocks);
    nucleation_parameters.validate(field_attributes, solve_blocks);
  }

  SpatialDiscretization<dim> spatial_discretization;
  TemporalDiscretization     temporal_discretization;
  BoundaryParameters<dim>    boundary_parameters;

  LinearSolveParameters    linear_solve_parameters;
  NonlinearSolveParameters nonlinear_solve_parameters;

  FieldOutputParameters   output_parameters;
  RestartOutputParameters restart_parameters;
  FieldInputParameters    input_parameters;

  MiscellaneousParameters misc_parameters;

  NucleationParameters nucleation_parameters;

  // TODO: This one needs to be fixed, but I don't want to touch it with a 9 foot pole.
  UserConstants<dim> user_constants;
};

PRISMS_PF_END_NAMESPACE
