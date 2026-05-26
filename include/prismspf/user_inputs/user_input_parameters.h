// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_block.h>

#include <prismspf/user_inputs/checkpoint_parameters.h>
#include <prismspf/user_inputs/constraint_parameters.h>
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
   * @brief Number of subsections to declare for certain field and solver parameters.
   */
  constexpr static unsigned int default_max_criteria = 5;

  /**
   * @brief Default Constructor.
   */
  UserInputParameters() = default;

  /**
   * @brief Constructor. Reads in user input parameters from file and loads them into
   * member variables.
   */
  explicit UserInputParameters(const std::string &file_name,
                               unsigned int       max_criteria = default_max_criteria);

  /**
   * @brief Tell parameter handler to expect the parameters by declaring them.
   */
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler,
                     unsigned int              max_criteria = default_max_criteria) const;

  /**
   * @brief Read the parameters from the parameter handler and assign them to the
   * appropriate member variables.
   */
  void
  assign_parameters(dealii::ParameterHandler &parameter_handler,
                    unsigned int              max_criteria = default_max_criteria);

  /**
   * @brief Ensure that the parameters are compatible with a set of fields and solvers.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks)
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

template <unsigned int dim>
inline UserInputParameters<dim>::UserInputParameters(const std::string &file_name,
                                                     unsigned int       max_criteria)
{
  // user_constants is a special little princess that needs to know the file contents
  // before we declare the parameters.
  user_constants.file_name = file_name;
  dealii::ParameterHandler parameter_handler;
  declare_parameters(parameter_handler, max_criteria);
  parameter_handler.parse_input(file_name);
  assign_parameters(parameter_handler, max_criteria);
}

template <unsigned int dim>
inline void
UserInputParameters<dim>::declare_parameters(dealii::ParameterHandler &parameter_handler,
                                             unsigned int              max_criteria) const
{
  // Perform and postprocessing of user inputs and run checks
  spatial_discretization.declare_parameters(parameter_handler, max_criteria);
  temporal_discretization.declare_parameters(parameter_handler);
  linear_solve_parameters.declare_parameters(parameter_handler, max_criteria);
  nonlinear_solve_parameters.declare_parameters(parameter_handler, max_criteria);
  output_parameters.declare_parameters(parameter_handler);
  checkpoint_parameters.declare_parameters(parameter_handler);
  boundary_parameters.declare_parameters(parameter_handler, max_criteria);
  nucleation_parameters.declare_parameters(parameter_handler);
  misc_parameters.declare_parameters(parameter_handler);
  load_ic_parameters.declare_parameters(parameter_handler);
  user_constants.declare_parameters(parameter_handler);
}

template <unsigned int dim>
inline void
UserInputParameters<dim>::assign_parameters(dealii::ParameterHandler &parameter_handler,
                                            unsigned int              max_criteria)
{
  // Perform and postprocessing of user inputs and run checks
  spatial_discretization.assign_parameters(parameter_handler, max_criteria);
  temporal_discretization.assign_parameters(parameter_handler);
  linear_solve_parameters.assign_parameters(parameter_handler, max_criteria);
  nonlinear_solve_parameters.assign_parameters(parameter_handler, max_criteria);
  output_parameters.assign_parameters(parameter_handler, temporal_discretization);
  checkpoint_parameters.assign_parameters(parameter_handler, temporal_discretization);
  boundary_parameters.assign_parameters(parameter_handler, max_criteria);
  nucleation_parameters.assign_parameters(parameter_handler);
  misc_parameters.assign_parameters(parameter_handler);
  load_ic_parameters.assign_parameters(parameter_handler);
  user_constants.assign_parameters(parameter_handler);
}

PRISMS_PF_END_NAMESPACE
