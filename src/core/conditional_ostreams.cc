#include <core/conditional_ostreams.h>

// NOLINTBEGIN

std::ofstream conditionalOStreams::summary_log_file("summary.log",
                                                    std::ios::out | std::ios::trunc);

dealii::ConditionalOStream conditionalOStreams::pout_summary_instance(
  summary_log_file,
  dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

const dealii::ConditionalOStream conditionalOStreams::pout_base(
  std::cout,
  dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

const dealii::ConditionalOStream conditionalOStreams::pout_verbose(
  std::cout,
#ifndef DEBUG
  false &&
#endif
    dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

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

dealii::ConditionalOStream &
conditionalOStreams::pout_summary()
{
  if (!summary_log_file.is_open())
    {
      throw std::runtime_error("summary.log file is not open.");
    }
  return pout_summary_instance;
}

// NOLINTEND