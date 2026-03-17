// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_group.h>

#include <prismspf/user_inputs/checkpoint_parameters.h>
#include <prismspf/user_inputs/constraint_parameters.h>
#include <prismspf/user_inputs/input_file_reader.h>
#include <prismspf/user_inputs/linear_solve_parameters.h>
#include <prismspf/user_inputs/load_initial_condition_parameters.h>
#include <prismspf/user_inputs/miscellaneous_parameters.h>
#include <prismspf/user_inputs/nonlinear_solve_parameters.h>
#include <prismspf/user_inputs/nucleation_parameters.h>
#include <prismspf/user_inputs/output_parameters.h>
#include <prismspf/user_inputs/spatial_discretization.h>
#include <prismspf/user_inputs/temporal_discretization.h>
#include <prismspf/user_inputs/user_constants.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
struct UserInputParameters
{
public:
  /**
   * @brief Default Constructor.
   */
  UserInputParameters() = default;
  /**
   * @brief Constructor. Reads in user input parameters from file and loads them into
   * member variables.
   */
  explicit UserInputParameters(const std::string &file_name);

  /**
   * @brief Ensure that the parameters are compatible with a set of fields and solvers.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveGroup>      &solve_groups)
  {
    // Perform and postprocessing of user inputs and run checks
    spatial_discretization.validate();
    temporal_discretization.validate();
    linear_solve_parameters.validate();
    nonlinear_solve_parameters.validate();
    output_parameters.validate();
    checkpoint_parameters.validate();
    boundary_parameters.validate();
    nucleation_parameters.validate();
    misc_parameters.postprocess_and_validate();
    load_ic_parameters.postprocess_and_validate();
  }

  /**
   * @brief Ensure that the parameters are compatible with a set of fields and solvers.
   */
  std::string
  parameter_summary()
  {
    // Print all the parameters to summary.log
    spatial_discretization.print_parameter_summary();
    temporal_discretization.print_parameter_summary();
    linear_solve_parameters.print_parameter_summary();
    nonlinear_solve_parameters.print_parameter_summary();
    output_parameters.print_parameter_summary();
    checkpoint_parameters.print_parameter_summary();
    boundary_parameters.print_parameter_summary();
    load_ic_parameters.print_parameter_summary();
    nucleation_parameters.print_parameter_summary();
    misc_parameters.print_parameter_summary();
    user_constants.print();
    return "";
  }

private:
  /**
   * @brief Assign the provided user inputs to parameters for anything related to the
   * spatial discretiziation.
   */
  void
  assign_spatial_discretization_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * @brief Assign the provided user inputs to parameters for anything related to the
   * temporal discretiziation.
   */
  void
  assign_temporal_discretization_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * @brief Assign the provided user inputs to parameters for anything related to linear
   * solves.
   */
  void
  assign_linear_solve_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * @brief Assign the provided user inputs to parameters for anything related to
   * nonlinear solves.
   */
  void
  assign_nonlinear_solve_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * @brief Assign the provided user inputs to parameters for anything related to
   * outputs.
   */
  void
  assign_output_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * @brief Assign the provided user inputs to parameters for anything related to
   * checkpoints.
   */
  void
  assign_checkpoint_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * @brief Assign the provided user inputs to parameters for anything related to
   * loading in initial condition.
   */
  void
  assign_load_initial_condition_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * @brief Assign the provided user inputs to parameters for anything related to
   * nucleation.
   */
  void
  assign_nucleation_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * @brief Assign the provided user inputs to parameters for anything related to
   * grain remapping and grain vtk load-in.
   */
  void
  assign_grain_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * @brief Assign the provided user inputs to parameters for anything related to
   * boundaries.
   */
  void
  assign_boundary_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * @brief Assign the provided user inputs to parameters for anything related to
   * miscellaneous parameters.
   */
  void
  assign_miscellaneous_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * @brief Assign the provided user constants.
   */
  void
  load_model_constants(const InputFileReader    &input_file_reader,
                       dealii::ParameterHandler &parameter_handler);

public:
  // Spatial discretization parameters
  SpatialDiscretization<dim> spatial_discretization;

  // Temporal discretization parameters
  TemporalDiscretization temporal_discretization;

  // Linear solve parameters
  LinearSolveParameters linear_solve_parameters;

  // Nonlinear solve parameters
  NonlinearSolveParameterSet nonlinear_solve_parameters;

  // Output parameters
  OutputParameters output_parameters;

  // Checkpoint parameters
  CheckpointParameters checkpoint_parameters;

  // Boundary parameters
  BoundaryParameters<dim> boundary_parameters;

  // Load IC parameters
  LoadInitialConditionParameters load_ic_parameters;

  // Nucleation parameters
  NucleationParameters nucleation_parameters;

  // Miscellaneous parameters
  MiscellaneousParameters misc_parameters;

  // User constants
  UserConstants<dim> user_constants;
};

PRISMS_PF_END_NAMESPACE
