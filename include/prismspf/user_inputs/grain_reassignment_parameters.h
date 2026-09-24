// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/exceptions.h>
#include <prismspf/core/simulation_timer.h>

#include <prismspf/user_inputs/parameter_base.h>
#include <prismspf/user_inputs/temporal_discretization.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include <climits>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that holds nucleation parameters.
 */
struct GrainReassignmentParameters : public ParameterBase
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  static void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int              n_subsections = Numbers::default_subsections);

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              n_subsections = Numbers::default_subsections) override;

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
  should_perform_grain_reassignment(unsigned int increment) const;

  // Whether to perform grain reassignment at all
  bool grain_reassignment_active = false;

  // The number of steps between grain reassignment
  unsigned int reassignment_period = UINT_MAX;

  // Whether to load the grain structure (WIP)
  bool load_grain_structure = false;

  // The distance between grains below which reassignment is required
  double exclusion_distance = 0.0;

  // The order parameter cutoff for identifying distinct grains
  double order_parameter_threshold = 0.01;
};

PRISMS_PF_END_NAMESPACE
