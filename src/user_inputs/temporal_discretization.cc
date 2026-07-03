#include <prismspf/core/solve_block.h>

#include <prismspf/user_inputs/parameter_base.h>
#include <prismspf/user_inputs/temporal_discretization.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

TemporalDiscretization::TemporalDiscretization(double _dt,
                                               double _initial_time,
                                               double _final_time)
  : dt(_dt)
  , initial_time(_initial_time)
  , final_time(_final_time)
  , n_increments(std::ceil(_final_time / _dt))
{}

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
  n_increments = (unsigned int) (parameter_handler.get_integer("final increment"));
  dt           = parameter_handler.get_double("time step");
  initial_time = parameter_handler.get_double("start time");
  final_time   = parameter_handler.get_double("end time");

  if (final_time > 0.0)
    {
      n_increments = std::ceil((final_time - initial_time) / dt);
    }
}

void
TemporalDiscretization::validate(
  [[maybe_unused]] const std::vector<FieldAttributes> &field_attributes,
  [[maybe_unused]] const std::vector<SolveBlock>      &solve_blocks) const
{
  double time_span = final_time - initial_time;

  AssertThrow(time_span > 0.0,
              dealii::ExcMessage("Final time must be greater than initial time."));
  AssertThrow(time_span / dt <= std::numeric_limits<unsigned int>::max(),
              dealii::ExcMessage("You seem to be taking more than 4 billion "
                                 "timesteps... That doesn't seem right. Right?"));
}

PRISMS_PF_END_NAMESPACE
