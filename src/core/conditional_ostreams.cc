// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/tee_stream.h>

#include <prismspf/config.h>

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <stdexcept>

PRISMS_PF_BEGIN_NAMESPACE

std::ofstream &
get_summary_log_file()
{
  static std::ofstream file("summary.log", std::ios::out | std::ios::trunc);
  if (!file.is_open())
    {
      throw std::runtime_error("Unable to open summary.log for writing.");
    }
  return file;
}

dealii::ConditionalOStream &
conditionalOStreams::pout_summary()
{
  static dealii::ConditionalOStream instance(get_summary_log_file(),
                                             dealii::Utilities::MPI::this_mpi_process(
                                               MPI_COMM_WORLD) == 0);
  return instance;
}

dealii::ConditionalOStream &
conditionalOStreams::pout_base()
{
  static TeeStream                  tee_stream(std::cout, get_summary_log_file());
  static dealii::ConditionalOStream instance(tee_stream,
                                             dealii::Utilities::MPI::this_mpi_process(
                                               MPI_COMM_WORLD) == 0);
  return instance;
}

dealii::ConditionalOStream &
conditionalOStreams::pout_verbose()
{
  static TeeStream                  tee_stream(std::cout, get_summary_log_file());
  static dealii::ConditionalOStream instance(tee_stream,
#ifndef DEBUG
                                             false &&
#endif
                                               dealii::Utilities::MPI::this_mpi_process(
                                                 MPI_COMM_WORLD) == 0);
  return instance;
}

PRISMS_PF_END_NAMESPACE
