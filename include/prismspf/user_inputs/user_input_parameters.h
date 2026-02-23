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

// TODO: convert to a simple struct

template <unsigned int dim>
class UserInputParameters
{
public:
  /**
   * @brief Constructor. Reads in user input parameters from file and loads them into
   * member variables.
   */
  UserInputParameters(InputFileReader          &input_file_reader,
                      dealii::ParameterHandler &parameter_handler);

  /**
   * @brief Return the spatial discretization parameters.
   */
  [[nodiscard]] const SpatialDiscretization<dim> &
  get_spatial_discretization() const
  {
    return spatial_discretization;
  }

  /**
   * @brief Return the spatial discretization parameters.
   */
  [[nodiscard]] SpatialDiscretization<dim> &
  get_spatial_discretization()
  {
    return spatial_discretization;
  }

  /**
   * @brief Return the temporal discretization parameters.
   */
  [[nodiscard]] const TemporalDiscretization &
  get_temporal_discretization() const
  {
    return temporal_discretization;
  }

  /**
   * @brief Return the temporal discretization parameters.
   */
  [[nodiscard]] TemporalDiscretization &
  get_temporal_discretization()
  {
    return temporal_discretization;
  }

  /**
   * @brief Return the linear solve parameters.
   */
  [[nodiscard]] const LinearSolveParameters &
  get_linear_solve_parameters() const
  {
    return linear_solve_parameters;
  }

  /**
   * @brief Return the linear solve parameters.
   */
  [[nodiscard]] LinearSolveParameters &
  get_linear_solve_parameters()
  {
    return linear_solve_parameters;
  }

  /**
   * @brief Return the nonlinear solve parameters.
   */
  [[nodiscard]] const NonlinearSolveParameterSet &
  get_nonlinear_solve_parameters() const
  {
    return nonlinear_solve_parameters;
  }

  /**
   * @brief Return the nonlinear solve parameters.
   */
  [[nodiscard]] NonlinearSolveParameterSet &
  get_nonlinear_solve_parameters()
  {
    return nonlinear_solve_parameters;
  }

  /**
   * @brief Return the output parameters.
   */
  [[nodiscard]] const OutputParameters &
  get_output_parameters() const
  {
    return output_parameters;
  }

  /**
   * @brief Return the output parameters.
   */
  [[nodiscard]] OutputParameters &
  get_output_parameters()
  {
    return output_parameters;
  }

  /**
   * @brief Return the checkpoint parameters.
   */
  [[nodiscard]] const CheckpointParameters &
  get_checkpoint_parameters() const
  {
    return checkpoint_parameters;
  }

  /**
   * @brief Return the checkpoint parameters.
   */
  [[nodiscard]] CheckpointParameters &
  get_checkpoint_parameters()
  {
    return checkpoint_parameters;
  }

  /**
   * @brief Return the boundary parameters.
   */
  [[nodiscard]] const BoundaryParameters<dim> &
  get_boundary_parameters() const
  {
    return boundary_parameters;
  }

  /**
   * @brief Return the boundary parameters.
   */
  [[nodiscard]] BoundaryParameters<dim> &
  get_boundary_parameters()
  {
    return boundary_parameters;
  }

  /**
   * @brief Return the load IC parameters.
   */
  [[nodiscard]] const LoadInitialConditionParameters &
  get_load_initial_condition_parameters() const
  {
    return load_ic_parameters;
  }

  /**
   * @brief Return the load IC parameters.
   */
  [[nodiscard]] LoadInitialConditionParameters &
  get_load_initial_condition_parameters()
  {
    return load_ic_parameters;
  }

  /**
   * @brief Return the nucleation parameters.
   */
  [[nodiscard]] const NucleationParameters &
  get_nucleation_parameters() const
  {
    return nucleation_parameters;
  }

  /**
   * @brief Return the nucleation parameters.
   */
  [[nodiscard]] NucleationParameters &
  get_nucleation_parameters()
  {
    return nucleation_parameters;
  }

  /**
   * @brief Return the miscellaneous parameters.
   */
  [[nodiscard]] const MiscellaneousParameters &
  get_miscellaneous_parameters() const
  {
    return misc_parameters;
  }

  /**
   * @brief Return the miscellaneous parameters.
   */
  [[nodiscard]] MiscellaneousParameters &
  get_miscellaneous_parameters()
  {
    return misc_parameters;
  }

  /**
   * @brief Return the user constants.
   */
  [[nodiscard]] const UserConstants<dim> &
  get_user_constants() const
  {
    return user_constants;
  }

  /**
   * @brief Return the user constants.
   */
  [[nodiscard]] UserConstants<dim> &
  get_user_constants()
  {
    return user_constants;
  }

  /**
   * @brief Ensure that the parameters are compatible with a set of fields and solvers.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveGroup>      &solve_groups)
  {
    // Perform and postprocessing of user inputs and run checks
    spatial_discretization.postprocess_and_validate();
    temporal_discretization.postprocess_and_validate();
    linear_solve_parameters.postprocess_and_validate();
    nonlinear_solve_parameters.postprocess_and_validate();
    output_parameters.postprocess_and_validate();
    checkpoint_parameters.postprocess_and_validate();
    boundary_parameters.postprocess_and_validate();
    nucleation_parameters.postprocess_and_validate();
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
