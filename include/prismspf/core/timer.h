// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/timer.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Timer class for PRISMS-PF
 *
 * This class keeps track of nested timing sections so we can print a structured summary
 * at the end of the simulation.
 *
 * Fortunately, `Caliper` handles this nicely (and more!) or us. `deal.II` doesn't do this
 * so we have to keep track of a few additional objects. The logical way to represent this
 * is with a tree node structure.
 */
class Timer
{
public:
  /**
   * @brief Constructor.
   */
  Timer();

  /**
   * @brief Start a new timer section.
   */
  static void
  start_section(const char *name);

  /**
   * @brief End the timer section.
   */
  static void
  end_section(const char *name);

  /**
   * @brief deal.II timer for the 0th MPI process
   */
  static dealii::TimerOutput &
  serial_timer();

  /**
   * @brief deal.II timer for parallel MPI process
   */
  static dealii::TimerOutput &
  parallel_timer();

  /**
   * @brief Print a sorted summary of the timed sections.
   */
  static void
  print_summary();
};

PRISMS_PF_END_NAMESPACE
