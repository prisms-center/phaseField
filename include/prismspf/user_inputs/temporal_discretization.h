// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/solve_block.h>

#include <prismspf/user_inputs/parameter_base.h>

#include <prismspf/config.h>

#include <cfloat>
#include <limits>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that holds temporal discretization parameters.
 */
struct TemporalDiscretization : public ParameterBase
{
  /**
   * @brief Constructor.
   */
  TemporalDiscretization() = default;

  /**
   * @brief Constructor.
   */
  TemporalDiscretization(double       _dt,
                         unsigned int _n_increments,
                         double       _initial_time = 0.0);

  /**
   * @brief Constructor.
   */
  TemporalDiscretization(double _dt, double _final_time, double _initial_time = 0.0);

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

  // Timestep
  double dt = 1.0;

  // Initial time
  double initial_time = 0.0;

  // Total number of increments
  unsigned int n_increments = 0;
};

PRISMS_PF_END_NAMESPACE
