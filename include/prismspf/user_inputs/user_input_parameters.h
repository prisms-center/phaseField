// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/user_inputs/boundary_parameters.h>
#include <prismspf/user_inputs/checkpoint_parameters.h>
#include <prismspf/user_inputs/linear_solve_parameters.h>
#include <prismspf/user_inputs/nonlinear_solve_parameters.h>
#include <prismspf/user_inputs/output_parameters.h>
#include <prismspf/user_inputs/spatial_discretization.h>
#include <prismspf/user_inputs/temporal_discretization.h>
#include <prismspf/user_inputs/user_constants.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

struct variableAttributes;

class inputFileReader;

template <unsigned int dim>
class userInputParameters
{
public:
  /**
   * \brief Constructor. Reads in user input parameters from file and loads them into
   * member variables.
   */
  userInputParameters(inputFileReader          &input_file_reader,
                      dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Return the variable attributes.
   */
  [[nodiscard]] const std::map<unsigned int, variableAttributes> &
  get_variable_attributes() const
  {
    return *var_attributes;
  }

  /**
   * \brief Return the spatial discretization parameters.
   */
  [[nodiscard]] const spatialDiscretization<dim> &
  get_spatial_discretization() const
  {
    return spatial_discretization;
  }

  /**
   * \brief Return the temporal discretization parameters.
   */
  [[nodiscard]] const temporalDiscretization &
  get_temporal_discretization() const
  {
    return temporal_discretization;
  }

  /**
   * \brief Return the linear solve parameters.
   */
  [[nodiscard]] const linearSolveParameters &
  get_linear_solve_parameters() const
  {
    return linear_solve_parameters;
  }

  /**
   * \brief Return the nonlinear solve parameters.
   */
  [[nodiscard]] const nonlinearSolveParameters &
  get_nonlinear_solve_parameters() const
  {
    return nonlinear_solve_parameters;
  }

  /**
   * \brief Return the output parameters.
   */
  [[nodiscard]] const outputParameters &
  get_output_parameters() const
  {
    return output_parameters;
  }

  /**
   * \brief Return the checkpoint parameters.
   */
  [[nodiscard]] const checkpointParameters &
  get_checkpoint_parameters() const
  {
    return checkpoint_parameters;
  }

  /**
   * \brief Return the boundary parameters.
   */
  [[nodiscard]] const boundaryParameters<dim> &
  get_boundary_parameters() const
  {
    return boundary_parameters;
  }

  /**
   * \brief Return the user constants.
   */
  [[nodiscard]] const userConstants<dim> &
  get_user_constants() const
  {
    return user_constants;
  }

private:
  /**
   * \brief Assign the provided user inputs to parameters for anything related to the
   * spatial discretiziation.
   */
  void
  assign_spatial_discretization_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to the
   * temporal discretiziation.
   */
  void
  assign_temporal_discretization_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to linear
   * solves.
   */
  void
  assign_linear_solve_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to
   * nonlinear solves.
   */
  void
  assign_nonlinear_solve_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to
   * outputs.
   */
  void
  assign_output_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to
   * checkpoints.
   */
  void
  assign_checkpoint_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to
   * loading in initial condition.
   */
  void
  assign_load_initial_condition_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to
   * nucleation.
   */
  void
  assign_nucleation_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to
   * grain remapping and grain vtk load-in.
   */
  void
  assign_grain_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to
   * boundaries.
   */
  void
  assign_boundary_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Assign the provided user constants.
   */
  void
  load_model_constants(const inputFileReader    &input_file_reader,
                       dealii::ParameterHandler &parameter_handler);

  // Variable attributes
  const std::map<unsigned int, variableAttributes> *var_attributes;

  // Spatial discretization parameters
  spatialDiscretization<dim> spatial_discretization;

  // Temporal discretization parameters
  temporalDiscretization temporal_discretization;

  // Linear solve paramters
  linearSolveParameters linear_solve_parameters;

  // Nonlinear solve parameters
  nonlinearSolveParameters nonlinear_solve_parameters;

  // Output parameters
  outputParameters output_parameters;

  // Checkpoint parameters
  checkpointParameters checkpoint_parameters;

  // Boundary parameters
  boundaryParameters<dim> boundary_parameters;

  // User constants
  userConstants<dim> user_constants;
};

PRISMS_PF_END_NAMESPACE
