// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/solve_group.h>

#include <prismspf/config.h>

#include <algorithm>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct that holds temporal discretization parameters.
 */
struct TemporalDiscretization
{
public:
  /**
   * @brief Construct from timestep and total number of increments,
   */
  explicit TemporalDiscretization(double _dt = 1.0, unsigned int _num_increments = 1)
    : dt(_dt)
    , num_increments(_num_increments)
  {}

  /**
   * @brief Construct from timestep and final time
   */
  TemporalDiscretization(double _dt, double final_time)
    : dt(_dt)
    , num_increments(uint((final_time + dt) / _dt))
  {}

  /**
   * @brief Postprocess and validate parameters.
   */
  void
  postprocess_and_validate();

  /**
   * @brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  /**
   * @brief Set the timestep.
   */
  void
  set_timestep(double _dt)
  {
    dt = _dt;
  }

  /**
   * @brief Get the timestep.
   */
  [[nodiscard]] double
  get_timestep() const
  {
    return dt;
  }

  /**
   * @brief Get the total number of increments.
   */
  [[nodiscard]] unsigned int
  get_num_increments() const
  {
    return num_increments;
  }

  /**
   * @brief Set the total number of increments.
   */
  void
  set_num_increments(unsigned int _total_increments)
  {
    num_increments = _total_increments;
  }

private:
  // Timestep
  double dt = 1.0;

  // Total number of increments
  unsigned int num_increments = 1;
};

inline void
TemporalDiscretization::print_parameter_summary() const
{
  ConditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Temporal Discretization\n"
    << "================================================\n"
    << "Timestep: " << dt << "\n"
    << "Total increments: " << num_increments << "\n\n"
    << std::flush;
}

PRISMS_PF_END_NAMESPACE
