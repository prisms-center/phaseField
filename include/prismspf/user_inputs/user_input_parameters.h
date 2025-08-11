// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef user_input_parameters_h
#define user_input_parameters_h

#include <prismspf/config.h>
#include <prismspf/user_inputs/boundary_parameters.h>
#include <prismspf/user_inputs/checkpoint_parameters.h>
#include <prismspf/user_inputs/input_file_reader.h>
#include <prismspf/user_inputs/linear_solve_parameters.h>
#include <prismspf/user_inputs/nonlinear_solve_parameters.h>
#include <prismspf/user_inputs/output_parameters.h>
#include <prismspf/user_inputs/spatial_discretization.h>
#include <prismspf/user_inputs/temporal_discretization.h>
#include <prismspf/user_inputs/user_constants.h>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
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
   * \brief Destructor.
   */
  ~userInputParameters() = default;

  // Variable attributes
  const std::map<unsigned int, variableAttributes> &var_attributes;

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
  load_model_constants(inputFileReader          &input_file_reader,
                       dealii::ParameterHandler &parameter_handler);
};

PRISMS_PF_END_NAMESPACE

#endif