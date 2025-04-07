// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/timer.h>

#include <prismspf/core/conditional_ostreams.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Timer class for PRISMS-PF
 */
class timer
{
public:
  /**
   * \brief Constructor.
   */
  timer();

  /**
   * \brief deal.II timer for the 0th MPI process
   */
  static dealii::TimerOutput &
  serial_timer();

  /**
   * \brief deal.II timer for parallel MPI process
   */
  static dealii::TimerOutput &
  parallel_timer();

  /**
   * \brief Print a sorted summary of the timed sections.
   */
  static void
  print_summary();
};

PRISMS_PF_END_NAMESPACE