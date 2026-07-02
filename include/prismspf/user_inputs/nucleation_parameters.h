// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/simulation_timer.h>

#include <prismspf/nucleation/nucleus.h>

#include <prismspf/user_inputs/parameter_base.h>
#include <prismspf/user_inputs/temporal_discretization.h>

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
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  predeclare(dealii::ParameterHandler &parameter_handler) const override;

  /**
   * @brief Assign the parameters from file.
   */
  void
  preassign(dealii::ParameterHandler &parameter_handler) override;

  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int max_criteria = Numbers::max_subsections) const override;

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              max_criteria = Numbers::max_subsections) override;

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override;

  /**
   * @brief Whether a given increment should attempt nucleation.
   */
  [[nodiscard]] bool
  should_attempt_nucleation(unsigned int increment) const;

  /**
   * @brief Check if a nucleus is still active based on its seed time and increment.
   * A nucleus is considered active if the current time or increment is within the greater
   * of the seeding time or seeding increments.
   * The exclusion distance and active refinement are only applied to active nuclei.
   */
  template <unsigned int dim>
  [[nodiscard]] bool
  check_active(const Nucleus<dim> &nucleus, const SimulationTimer &time_info) const;

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

PRISMS_PF_END_NAMESPACE
