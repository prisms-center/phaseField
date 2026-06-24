// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/solve_block.h>

#include <prismspf/config.h>

#include "prismspf/user_inputs/user_input_parameters.h"

#include <cfloat>
#include <limits>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that holds temporal discretization parameters.
 */
struct TemporalDiscretization : public ParameterBase
{
public:
  /**
   * @brief Constructor.
   */
  TemporalDiscretization() = default;

  /**
   * @brief Constructor.
   */
  TemporalDiscretization(double _dt, double _initial_time, double _final_time)
    : dt(_dt)
    , initial_time(_initial_time)
    , final_time(_final_time)
    , num_increments(std::ceil(_final_time / _dt))
  {}

  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int max_criteria = Numbers::max_subsections) const override
  {
    parameter_handler.declare_entry("final increment",
                                    "0",
                                    dealii::Patterns::Integer(0, INT_MAX),
                                    "The final increment for the simulation.");
    declare_aliases(parameter_handler,
                    "final increment",
                    {"number steps",
                     "number of steps",
                     "num steps",
                     "n steps",
                     "last increment",
                     "end increment",
                     "final_increment",
                     "number_steps",
                     "number_of_steps",
                     "num_steps",
                     "n_steps",
                     "last_increment",
                     "end_increment",
                     "steps",
                     "iterations",
                     "final iteration",
                     "max iteration",
                     "increments"});

    parameter_handler.declare_entry("time step",
                                    "0.0",
                                    dealii::Patterns::Double(0.0, DBL_MAX),
                                    "The time step size for the simulation.");
    declare_aliases(parameter_handler, "time step", {"timestep", "time_step", "dt"});

    parameter_handler.declare_entry("start time",
                                    "0.0",
                                    dealii::Patterns::Double(0.0, DBL_MAX),
                                    "The value of simulated time where the simulation "
                                    "begins.");
    declare_aliases(parameter_handler,
                    "start time",
                    {"start_time", "t_0", "t0", "begin time", "begin_time"});

    parameter_handler.declare_entry("end time",
                                    "0.0",
                                    dealii::Patterns::Double(0.0, DBL_MAX),
                                    "The value of simulated time where the simulation "
                                    "ends. Overrides final increment if greater than 0.");
    declare_aliases(parameter_handler,
                    "end time",
                    {"end_time", "t_f", "tf", "final time", "final_time"});
  };

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              max_criteria = Numbers::max_subsections) override
  {
    n_increments = (unsigned int) (parameter_handler.get_integer("final increment"));
    dt           = parameter_handler.get_double("time step");
    initial_time = parameter_handler.get_double("start time");
    final_time   = parameter_handler.get_double("end time");

    if (final_time > 0.0)
      {
        n_increments = std::ceil((final_time - initial_time) / dt);
      }
  };

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override
  {
    double time_span = final_time - initial_time;

    AssertThrow(time_span > 0.0,
                dealii::ExcMessage("Final time must be greater than initial time."));
    AssertThrow(time_span / dt <= std::numeric_limits<unsigned int>::max(),
                dealii::ExcMessage("You seem to be taking more than 4 billion "
                                   "timesteps... That doesn't seem right. Right?"));
  };

  // Timestep
  double dt = 1.0;

  // Initial time
  double initial_time = 0.0;

  // Final time
  double final_time = 1.0;

  // Total number of increments
  unsigned int n_increments = 1;
};

PRISMS_PF_END_NAMESPACE
