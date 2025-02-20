// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>

#include <ios>
#include <iostream>
#include <mpi.h>
#include <stdexcept>

PRISMS_PF_BEGIN_NAMESPACE

// NOLINTBEGIN

std::ofstream conditionalOStreams::summary_log_file("summary.log",
                                                    std::ios::out | std::ios::trunc);

// NOLINTEND

dealii::ConditionalOStream &
conditionalOStreams::pout_summary()
{
  static dealii::ConditionalOStream instance(summary_log_file,
                                             dealii::Utilities::MPI::this_mpi_process(
                                               MPI_COMM_WORLD) == 0);
  return instance;
}

dealii::ConditionalOStream &
conditionalOStreams::pout_base()
{
  static dealii::ConditionalOStream instance(std::cout,
                                             dealii::Utilities::MPI::this_mpi_process(
                                               MPI_COMM_WORLD) == 0);
  return instance;
}

dealii::ConditionalOStream &
conditionalOStreams::pout_verbose()
{
  static dealii::ConditionalOStream instance(std::cout,
#ifndef DEBUG
                                             false &&
#endif
                                               dealii::Utilities::MPI::this_mpi_process(
                                                 MPI_COMM_WORLD) == 0);
  return instance;
}

conditionalOStreams::conditionalOStreams()
{
  if (!summary_log_file.is_open())
    {
      throw std::runtime_error("Unable to open summary.log for writing.");
    }
}

conditionalOStreams::~conditionalOStreams()
{
  if (summary_log_file.is_open())
    {
      summary_log_file.close();
    }
}

PRISMS_PF_END_NAMESPACE