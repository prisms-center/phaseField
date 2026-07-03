#include <prismspf/core/solve_block.h>

#include <prismspf/user_inputs/parameter_base.h>
#include <prismspf/user_inputs/temporal_discretization.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

TemporalDiscretization::TemporalDiscretization(double       _dt,
                                               unsigned int _n_increments,
                                               double       _initial_time)
  : dt(_dt)
  , initial_time(_initial_time)
  , n_increments(_n_increments)
{}

TemporalDiscretization::TemporalDiscretization(double _dt,
                                               double _final_time,
                                               double _initial_time)
  : dt(_dt)
  , initial_time(_initial_time)
  , n_increments(_final_time == _initial_time
                   ? 0
                   : (unsigned int) std::ceil((_final_time - _initial_time) / _dt))
{
  AssertThrow(initial_time <= _final_time,
              dealii::ExcMessage(
                "Initial time must be less than or equal to final time."));
}

void
TemporalDiscretization::declare(dealii::ParameterHandler &parameter_handler,
                                unsigned int              n_subsections)
{
  parameter_handler.declare_entry("final increment",
                                  "0",
                                  dealii::Patterns::Integer(0, INT_MAX),
                                  "The final increment for the simulation.");
  declare_aliases(parameter_handler,
                  "final increment",
                  std::vector {"number steps",
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
  declare_aliases(parameter_handler,
                  "time step",
                  std::vector {"timestep", "time_step", "dt"});

  parameter_handler.declare_entry("start time",
                                  "0.0",
                                  dealii::Patterns::Double(0.0, DBL_MAX),
                                  "The value of simulated time where the simulation "
                                  "begins.");
  declare_aliases(parameter_handler,
                  "start time",
                  std::vector {"start_time", "t_0", "t0", "begin time", "begin_time"});

  parameter_handler.declare_entry("end time",
                                  "0.0",
                                  dealii::Patterns::Double(0.0, DBL_MAX),
                                  "The value of simulated time where the simulation "
                                  "ends. Overrides final increment if greater than 0.");
  declare_aliases(parameter_handler,
                  "end time",
                  std::vector {"end_time", "t_f", "tf", "final time", "final_time"});
}

void
TemporalDiscretization::assign(dealii::ParameterHandler &parameter_handler,
                               unsigned int              n_subsections)
{
  n_increments      = (unsigned int) parameter_handler.get_integer("final increment");
  dt                = parameter_handler.get_double("time step");
  initial_time      = parameter_handler.get_double("start time");
  double final_time = parameter_handler.get_double("end time");

  if (final_time > 0.0)
    {
      AssertThrow(initial_time <= final_time,
                  dealii::ExcMessage(
                    "Initial time must be less than or equal to final time."));
      n_increments = std::ceil((final_time - initial_time) / dt);
    }
}

void
TemporalDiscretization::validate(
  [[maybe_unused]] const std::vector<FieldAttributes> &field_attributes,
  [[maybe_unused]] const std::vector<SolveBlock>      &solve_blocks) const
{
  AssertThrow(n_increments == 0 || dt > 0.0,
              dealii::ExcMessage(
                "Time step must be greater than 0 for transient problems."));
}

PRISMS_PF_END_NAMESPACE
