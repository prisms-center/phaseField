// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>

#include <prismspf/user_inputs/temporal_discretization.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>
#include <prismspf/nucleation/nucleus.h>

#include "prismspf/core/variable_attributes.h"

#include <climits>
#include <set>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that holds nucleation parameters.
 */
struct NucleationParameters
{
public:
  /**
   * @brief Return if the increment should be nucleationted.
   */
  [[nodiscard]] bool
  should_attempt_nucleation(unsigned int increment) const;

  /**
   * @brief Postprocess and validate parameters.
   */
  void
  postprocess_and_validate(
    const std::map<unsigned int, VariableAttributes> &var_attributes);

  /**
   * @brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  /**
   * @brief Set the print nucleation period
   */
  void
  set_nucleation_period(const unsigned int &_nucleation_period)
  {
    nucleation_period = _nucleation_period;
  }

  /**
   * @brief Get the nucleation period
   */
  [[nodiscard]] unsigned int
  get_nucleation_period() const
  {
    return nucleation_period;
  }

  /**
   * @brief Set the nucleus exclusion distance
   */
  void
  set_exclusion_distance(const double &_nucleus_exclusion_distance)
  {
    AssertThrow(_nucleus_exclusion_distance >= 0.0,
                dealii::ExcMessage("Nucleus exclusion distance must be non-negative."));
    nucleus_exclusion_distance = _nucleus_exclusion_distance;
  }

  /**
   * @brief Get the nucleus exclusion distance
   */
  [[nodiscard]] double
  get_exclusion_distance() const
  {
    return nucleus_exclusion_distance;
  }

  /**
   * @brief Set the nucleus exclusion distance
   */
  void
  set_same_field_exclusion_distance(const double &exclusion_distance)
  {
    AssertThrow(exclusion_distance >= 0.0,
                dealii::ExcMessage("Nucleus exclusion distance must be non-negative."));
    same_field_nucleus_exclusion_distance = exclusion_distance;
  }

  /**
   * @brief Get the nucleus exclusion distance
   */
  [[nodiscard]] double
  get_same_field_exclusion_distance() const
  {
    return same_field_nucleus_exclusion_distance;
  }

  /**
   * @brief Set the refinement radius
   */
  void
  set_refinement_radius(const double &_refinement_radius)
  {
    AssertThrow(_refinement_radius >= 0.0,
                dealii::ExcMessage("Nucleation refinement radius must be non-negative."));
    refinement_radius = _refinement_radius;
  }

  /**
   * @brief Get the refinement radius
   */
  [[nodiscard]] double
  get_refinement_radius() const
  {
    return refinement_radius;
  }

  /**
   * @brief Set the seeding time
   */
  void
  set_seeding_time(const double &_seeding_time)
  {
    seeding_time = _seeding_time;
  }

  /**
   * @brief Get the seeding time
   */
  [[nodiscard]] double
  get_seeding_time() const
  {
    return seeding_time;
  }

  /**
   * @brief Set the seeding increments
   */
  void
  set_seeding_increments(const unsigned int &_seeding_increments)
  {
    seeding_increments = _seeding_increments;
  }

  /**
   * @brief Get the seeding increments
   */
  [[nodiscard]] unsigned int
  get_seeding_increments() const
  {
    return seeding_increments;
  }

  /**
   * @brief Check if a nucleus is still active based on its seed time and increment.
   * A nucleus is considered active if the current time or increment is within the greater
   * of the seeding time or seeding increments.
   * The exclusion distance and active refinement are only applied to active nuclei.
   */
  template <unsigned int dim>
  [[nodiscard]] bool
  check_active(const Nucleus<dim> &nucleus, const TemporalDiscretization &time_info) const
  {
    return nucleus.seed_increment <= time_info.get_increment() &&
           ((time_info.get_increment() - nucleus.seed_increment) < seeding_increments ||
            time_info.get_time() - nucleus.seed_time < seeding_time);
  }

  /**
   * @brief Whether to print timing information with nucleation
   */
  void
  set_print_timing_with_nucleation(const bool &_print_timing_with_nucleation)
  {
    print_timing_with_nucleation = _print_timing_with_nucleation;
  }

  /**
   * @brief Whether a postprocessed nucleation rate exists
   */
  [[nodiscard]] bool
  postprocessed_nucleation_rate_exists() const
  {
    return pp_nucleation_rate_exists;
  }

private:
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

  // Whether to print timing information with nucleation
  // TODO (landinjm): Implement this.
  bool print_timing_with_nucleation = false;
};

inline bool
NucleationParameters::should_attempt_nucleation(unsigned int increment) const
{
  return !bool(increment % nucleation_period);
}

inline void
NucleationParameters::postprocess_and_validate(
  const std::map<unsigned int, VariableAttributes> &var_attributes)
{
  // Check if the nucleation period is valid
  AssertThrow(nucleation_period > 0,
              dealii::ExcMessage("Nucleation period must be positive."));

  // Check if the postprocessed nucleation rate exists
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.is_postprocess() && variable.is_nucleation_rate())
        {
          pp_nucleation_rate_exists = true;
          break;
        }
    }
}

inline void
NucleationParameters::print_parameter_summary() const
{
  ConditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Nucleation Parameters\n"
    << "================================================\n"
    << "Nucleation period: " << nucleation_period << "\n"
    << "Print nucleation timing info: " << bool_to_string(print_timing_with_nucleation)
    << "\n\n"
    << std::flush;
}

PRISMS_PF_END_NAMESPACE
