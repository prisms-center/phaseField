// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/user_inputs/nucleation_parameters.h>
#include <prismspf/user_inputs/parameter_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
NucleationParameters::declare(dealii::ParameterHandler &parameter_handler,
                              unsigned int              n_subsections)
{
  parameter_handler.enter_subsection("nucleation");
  {
    declare_entry(parameter_handler,
                  "nucleus exclusion distance",
                  "0.0",
                  Patterns::PositiveDouble(),
                  "The minimum distance between seeded nuclei for all nucleating fields.",
                  {"exclusion distance", "nucleus exclusion radius", "exclusion radius"});
    declare_entry(
      parameter_handler,
      "same field nucleus exclusion distance",
      "0.0",
      Patterns::PositiveDouble(),
      "The minimum distance between seeded nuclei within the same nucleating field.",
      {"same field exclusion distance",
       "same field nucleus exclusion radius",
       "same field exclusion radius"});
    declare_entry(parameter_handler,
                  "nucleation stride",
                  "2147483647",
                  Patterns::PositiveInteger(),
                  "The number of steps between nucleation attempts.",
                  {"nucleation period",
                   "steps between nucleation",
                   "number of steps between nucleation",
                   "n steps between nucleation",
                   "num steps between nucleation"});
    declare_entry(parameter_handler,
                  "refinement radius",
                  "0.0",
                  Patterns::PositiveDouble(),
                  "The radius around a nucleus in which AMR is applied.");
    declare_entry(
      parameter_handler,
      "seeding time",
      "0.0",
      Patterns::PositiveDouble(),
      "The time duration over which nuclei are considered \"active\" and refinement and "
      "exclusion zones are applied. Same as \"seeding increments\" but in time.",
      {"t_seed"});
    declare_entry(
      parameter_handler,
      "seeding increments",
      "1",
      Patterns::PositiveInteger(),
      "The number of increments over which nuclei are considered \"active\" and "
      "refinement and exclusion zones are applied. Same as \"seeding time\" but in "
      "increments.",
      {"seeding steps", "seeding iterations"});
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

    nucleation_period = (unsigned int) parameter_handler.get_integer("nucleation stride");

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
