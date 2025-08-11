// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef conditional_ostreams_h
#define conditional_ostreams_h

#include <deal.II/base/conditional_ostream.h>

#include <prismspf/config.h>

#include <fstream>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief A class that allows printing to different output streams that are classified
 * based on their verbosity. For now, this consists of two stream the release and debug.
 * The debug stream provides more information that may be useful when debugging.
 */
class conditionalOStreams
{
public:
  /**
   * \brief Constructor.
   */
  conditionalOStreams();

  /**
   * \brief Destructor.
   */
  ~conditionalOStreams();

  /**
   * \brief Generic parallel output stream. Used for essential information in release and
   * debug mode.
   */
  static dealii::ConditionalOStream &
  pout_base();

  /**
   * \brief Verbose parallel output stream. Used for additional information in debug mode.
   */
  static dealii::ConditionalOStream &
  pout_verbose();

  /**
   * \brief Log output stream for writing a summary.log file.
   */
  static dealii::ConditionalOStream &
  pout_summary();

private:
  // summary.log file
  static std::ofstream summary_log_file;
};

PRISMS_PF_END_NAMESPACE

#endif