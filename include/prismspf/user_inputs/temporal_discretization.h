// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/solve_block.h>

#include <prismspf/config.h>

#include <cfloat>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that holds temporal discretization parameters.
 */
struct TemporalDiscretization
{
public:
  /**
   * @brief Construct from timestep and total number of increments,
   */
  explicit TemporalDiscretization(double _dt = 1.0, unsigned int _num_increments = 1)
    : dt(_dt)
    , num_increments(_num_increments)
  {}

  /**
   * @brief Construct from timestep and final time
   */
  TemporalDiscretization(double _dt, double final_time)
    : dt(_dt)
    , num_increments(static_cast<unsigned int>((final_time + dt) / _dt))
  {}

  /**
   * @brief Postprocess and validate parameters.
   */
  void
  validate()
  {}

  /**
   * @brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler) const;

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  void
  assign_parameters(dealii::ParameterHandler &parameter_handler);

  // Timestep
  double dt = 1.0;

  // Total number of increments
  unsigned int num_increments = 0;
};

inline void
TemporalDiscretization::print_parameter_summary() const
{
  ConditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Temporal Discretization\n"
    << "================================================\n"
    << "Timestep: " << dt << "\n"
    << "Total increments: " << num_increments << "\n\n"
    << std::flush;
}

inline void
TemporalDiscretization::declare_parameters(
  dealii::ParameterHandler &parameter_handler) const
{
  parameter_handler.declare_entry("final increment",
                                  "0",
                                  dealii::Patterns::Integer(0, INT_MAX),
                                  "The final increment for the simulation.");
  parameter_handler.declare_entry("time step",
                                  "0.0",
                                  dealii::Patterns::Double(0.0, DBL_MAX),
                                  "The time step size for the simulation.");
  parameter_handler.declare_entry("end time",
                                  "0.0",
                                  dealii::Patterns::Double(0.0, DBL_MAX),
                                  "The value of simulated time where the simulation "
                                  "ends. Overrides final increment if greater than 0.");

  parameter_handler.declare_alias("final increment", "number steps");
  parameter_handler.declare_alias("final increment", "number of steps");
  parameter_handler.declare_alias("final increment", "num steps");
  parameter_handler.declare_alias("final increment", "n steps");
  parameter_handler.declare_alias("final increment", "last increment");
  parameter_handler.declare_alias("final increment", "end increment");
  parameter_handler.declare_alias("final increment", "final_increment");
  parameter_handler.declare_alias("final increment", "number_steps");
  parameter_handler.declare_alias("final increment", "number_of_steps");
  parameter_handler.declare_alias("final increment", "num_steps");
  parameter_handler.declare_alias("final increment", "n_steps");
  parameter_handler.declare_alias("final increment", "last_increment");
  parameter_handler.declare_alias("final increment", "end_increment");
  parameter_handler.declare_alias("final increment", "steps"); //
  parameter_handler.declare_alias("final increment", "iterations");
  parameter_handler.declare_alias("final increment", "increments");

  parameter_handler.declare_alias("time step", "timestep");
  parameter_handler.declare_alias("time step", "time_step");
  parameter_handler.declare_alias("time step", "dt");

  parameter_handler.declare_alias("end time", "final time");
  parameter_handler.declare_alias("end time", "end_time");
}

inline void
TemporalDiscretization::assign_parameters(dealii::ParameterHandler &parameter_handler)
{
  dt = parameter_handler.get_double("time step");
  num_increments =
    static_cast<unsigned int>(parameter_handler.get_integer("final increment"));
  double final_time = parameter_handler.get_double("end time");
  if (final_time > 0.0)
    {
      num_increments = std::ceil(final_time / dt);
    }
}

PRISMS_PF_END_NAMESPACE
