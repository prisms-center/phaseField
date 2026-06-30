// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/simulation_timer.h>

#include <prismspf/nucleation/nucleus.h>

#include <prismspf/user_inputs/temporal_discretization.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include <climits>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that holds nucleation parameters.
 */
struct NucleationParameters : public ParameterBase
{
public:
  /**
   * @brief Whether a given increment should attempt nucleation.
   */
  [[nodiscard]] bool
  should_attempt_nucleation(unsigned int increment) const
  {
    return increment % nucleation_period == 0;
  };

  /**
   * @brief Check if a nucleus is still active based on its seed time and increment.
   * A nucleus is considered active if the current time or increment is within the greater
   * of the seeding time or seeding increments.
   * The exclusion distance and active refinement are only applied to active nuclei.
   */
  template <unsigned int dim>
  [[nodiscard]] bool
  check_active(const Nucleus<dim> &nucleus, const SimulationTimer &time_info) const
  {
    return nucleus.seed_increment <= time_info.get_increment() &&
           ((time_info.get_increment() - nucleus.seed_increment) < seeding_increments ||
            time_info.get_time() - nucleus.seed_time < seeding_time);
  }

  // The number of steps between nucleationting relevant information to screen
  unsigned int nucleation_period = UINT_MAX;

  // Whether a postprocess field is being used for nucleation
  bool pp_nucleation_rate_exists = false;

  // The radius around a nucleus to exclude other nuclei
  double nucleus_exclusion_distance = 0.0;

  // The radius around a nucleus to exclude other nuclei in the same field
  double same_field_nucleus_exclusion_distance = 0.0;

  // The radius around a nucleus to refine the mesh
  double refinement_radius = 0.0;

  // Seeding time
  double seeding_time = 0.0;

  // Seeding increments
  unsigned int seeding_increments = 1;
};

inline void
NucleationParameters::declare_parameters(
  dealii::ParameterHandler &parameter_handler) const
{
  parameter_handler.enter_subsection("nucleation");
  {
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
    { // Declare aliases for the parameters
      //============================================================================================
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
      //
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
  }
  parameter_handler.leave_subsection();
}

inline void
NucleationParameters::assign_parameters(dealii::ParameterHandler &parameter_handler)
{
  parameter_handler.enter_subsection("nucleation");
  {
    set_exclusion_distance(parameter_handler.get_double("nucleus exclusion distance"));

    set_same_field_exclusion_distance(
      parameter_handler.get_double("same field nucleus exclusion distance"));

    set_nucleation_period(
      static_cast<unsigned int>(parameter_handler.get_integer("nucleation period")));

    set_refinement_radius(parameter_handler.get_double("refinement radius"));

    set_seeding_time(parameter_handler.get_double("seeding time"));

    set_seeding_increments(
      static_cast<unsigned int>(parameter_handler.get_integer("seeding increments")));
  }
  parameter_handler.leave_subsection();
}

PRISMS_PF_END_NAMESPACE
