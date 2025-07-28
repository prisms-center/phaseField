// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/conditional_ostream.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief A class that allows printing to different output streams that are classified
 * based on their verbosity. For now, this consists of two stream the release and debug.
 * The debug stream provides more information that may be useful when debugging.
 */
class ConditionalOStreams
{
public:
  /**
   * @brief Constructor.
   */
  ConditionalOStreams() = default;

  /**
   * @brief Destructor.
   */
  ~ConditionalOStreams() = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so ostream instances aren't copied.
   */
  ConditionalOStreams(const ConditionalOStreams &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so ostream instances aren't copied.
   */
  ConditionalOStreams &
  operator=(const ConditionalOStreams &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so ostream instances aren't moved.
   */
  ConditionalOStreams(ConditionalOStreams &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so ostream instances aren't moved.
   */
  ConditionalOStreams &
  operator=(ConditionalOStreams &&solver_base) noexcept = delete;

  /**
   * @brief Generic parallel output stream. Used for essential information in release and
   * debug mode.
   */
  static dealii::ConditionalOStream &
  pout_base();

  /**
   * @brief Verbose parallel output stream. Used for additional information in debug mode.
   */
  static dealii::ConditionalOStream &
  pout_verbose();

  /**
   * @brief Log output stream for writing a summary.log file.
   */
  static dealii::ConditionalOStream &
  pout_summary();
};

PRISMS_PF_END_NAMESPACE