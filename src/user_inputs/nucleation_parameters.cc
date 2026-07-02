#include <prismspf/user_inputs/nucleation_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
NucleationParameters::declare(dealii::ParameterHandler &parameter_handler,
                              unsigned int              n_subsections) const
{
  parameter_handler.enter_subsection("nucleation");
  {
    // TODO: See bounds for these
    parameter_handler.declare_entry("nucleus exclusion distance",
                                    "0.0",
                                    dealii::Patterns::Double(),
                                    "The minimum distance between nuclei.");
    parameter_handler.declare_entry("same field nucleus exclusion distance",
                                    "0.0",
                                    dealii::Patterns::Double(),
                                    "The minimum distance between nuclei.");
    parameter_handler.declare_entry(
      "nucleation period",
      "2147483647",
      dealii::Patterns::Integer(1),
      "The number of increments between nucleation attempts.");
    parameter_handler.declare_entry(
      "refinement radius",
      "0.0",
      dealii::Patterns::Double(0.0),
      "The radius around a nucleus in which AMR is applied.");
    parameter_handler.declare_entry(
      "seeding time",
      "0.0",
      dealii::Patterns::Double(0.0),
      "The time duration over which nuclei are considered \"active\" and refinement and "
      "exclusion zones are applied. Same as \"seeding increments\" but in time.");
    parameter_handler.declare_entry(
      "seeding increments",
      "1",
      dealii::Patterns::Integer(1, INT_MAX),
      "The number of increments over which nuclei are considered \"active\" and "
      "refinement and exclusion zones are applied. Same as \"seeding time\" but in "
      "increments.");

    // TODO: Use helper functions for these
    parameter_handler.declare_alias("nucleus exclusion distance",
                                    "nucleus_exclusion_distance");
    parameter_handler.declare_alias("nucleus exclusion distance",
                                    "nucleus exclusion radius");
    parameter_handler.declare_alias("nucleus exclusion distance",
                                    "nucleus_exclusion_radius");
    parameter_handler.declare_alias("nucleus exclusion distance", "exclusion distance");
    parameter_handler.declare_alias("nucleus exclusion distance", "exclusion_distance");
    parameter_handler.declare_alias("nucleus exclusion distance", "exclusion radius");
    parameter_handler.declare_alias("nucleus exclusion distance", "exclusion_radius");
    parameter_handler.declare_alias("same field nucleus exclusion distance",
                                    "same_field_nucleus_exclusion_distance");
    parameter_handler.declare_alias("same field nucleus exclusion distance",
                                    "same field nucleus exclusion radius");
    parameter_handler.declare_alias("same field nucleus exclusion distance",
                                    "same_field_nucleus_exclusion_radius");
    parameter_handler.declare_alias("same field nucleus exclusion distance",
                                    "same field exclusion distance");
    parameter_handler.declare_alias("same field nucleus exclusion distance",
                                    "same_field_exclusion_distance");
    parameter_handler.declare_alias("same field nucleus exclusion distance",
                                    "same field exclusion radius");
    parameter_handler.declare_alias("same field nucleus exclusion distance",
                                    "same_field_exclusion_radius");
  }
  parameter_handler.leave_subsection();
}

void
NucleationParameters::assign(dealii::ParameterHandler &parameter_handler,
                             unsigned int              n_subsections)
{
  parameter_handler.enter_subsection("nucleation");
  {
    nucleus_exclusion_distance =
      parameter_handler.get_double("nucleus exclusion distance");

    same_field_nucleus_exclusion_distance =
      parameter_handler.get_double("same field nucleus exclusion distance");

    nucleation_period = (unsigned int) parameter_handler.get_integer("nucleation period");

    refinement_radius = parameter_handler.get_double("refinement radius");

    seeding_time = parameter_handler.get_double("seeding time");

    seeding_increments =
      (unsigned int) parameter_handler.get_integer("seeding increments");
  }
  parameter_handler.leave_subsection();
}

void
NucleationParameters::validate(const std::vector<FieldAttributes> &field_attributes,
                               const std::vector<SolveBlock>      &solve_blocks) const
{
  // TODO: Do this later
}

bool
NucleationParameters::should_attempt_nucleation(unsigned int increment) const
{
  return increment % nucleation_period == 0;
}

template <unsigned int dim>
bool
NucleationParameters::check_active(const Nucleus<dim>    &nucleus,
                                   const SimulationTimer &time_info) const
{
  return nucleus.seed_increment <= time_info.get_increment() &&
         ((time_info.get_increment() - nucleus.seed_increment) < seeding_increments ||
          time_info.get_time() - nucleus.seed_time < seeding_time);
}

#include "user_inputs/nucleation_parameters.inst"

PRISMS_PF_END_NAMESPACE
