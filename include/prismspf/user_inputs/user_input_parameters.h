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

#include <prismspf/utilities/assert.h>

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
                               unsigned int n_subsections = Numbers::default_subsections)
  {
    // user_constants is a special little princess that needs to know the file contents
    // before we declare the parameters.
    user_constants.file_name = file_name;

    DEBUG_ASSERT(1 == 2, "Debug assert");
    ASSERT(1 == 4, "assert");

    dealii::ParameterHandler parameter_handler;

    declare(parameter_handler, n_subsections);
    parameter_handler.parse_input(file_name);
    assign(parameter_handler, n_subsections);
  };

  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int              n_subsections = Numbers::default_subsections) const
  {
    SpatialDiscretization<dim>::declare(parameter_handler, n_subsections);
    TemporalDiscretization::declare(parameter_handler, n_subsections);
    BoundaryParameters::declare(parameter_handler, n_subsections);

    LinearSolveParameters::declare(parameter_handler, n_subsections);
    NonlinearSolveParameters::declare(parameter_handler, n_subsections);

    FieldOutputParameters::declare(parameter_handler, n_subsections);
    RestartOutputParameters::declare(parameter_handler, n_subsections);
    FieldInputParameters::declare(parameter_handler, n_subsections);

    MiscellaneousParameters::declare(parameter_handler, n_subsections);
    NucleationParameters::declare(parameter_handler, n_subsections);

    user_constants.declare_parameters(parameter_handler);
  };

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              n_subsections = Numbers::default_subsections) override
  {
    spatial_discretization.assign(parameter_handler, n_subsections);
    temporal_discretization.assign(parameter_handler, n_subsections);
    boundary_parameters.assign<dim>(parameter_handler, n_subsections);

    linear_solve_parameters.assign(parameter_handler, n_subsections);
    nonlinear_solve_parameters.assign(parameter_handler, n_subsections);

    const auto n_increments = temporal_discretization.n_increments;

    output_parameters.assign(parameter_handler, n_increments, n_subsections);
    restart_parameters.assign(parameter_handler, n_increments, n_subsections);
    input_parameters.assign(parameter_handler, n_subsections);

    misc_parameters.assign(parameter_handler, n_subsections);
    nucleation_parameters.assign(parameter_handler, n_subsections);

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
  BoundaryParameters         boundary_parameters;

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
