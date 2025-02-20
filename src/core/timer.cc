// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/timer.h>

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/timer.h>

#include <iostream>
#include <mpi.h>

PRISMS_PF_BEGIN_NAMESPACE

dealii::TimerOutput &
timer::serial_timer()
{
  static dealii::TimerOutput instance(conditionalOStreams::pout_base(),
                                      dealii::TimerOutput::never,
                                      dealii::TimerOutput::wall_times);

  return instance;
}

dealii::TimerOutput &
timer::parallel_timer()
{
  static dealii::TimerOutput instance(MPI_COMM_WORLD,
                                      std::cout,
                                      dealii::TimerOutput::never,
                                      dealii::TimerOutput::wall_times);

  return instance;
}

void
timer::print_summary()
{
  // Get the timer output for the serial and parallel timers
  const auto serial_n_calls =
    serial_timer().get_summary_data(dealii::TimerOutput::OutputData::n_calls);
  const auto serial_wall_times =
    serial_timer().get_summary_data(dealii::TimerOutput::OutputData::total_wall_time);
  const auto parallel_n_calls =
    parallel_timer().get_summary_data(dealii::TimerOutput::OutputData::n_calls);
  const auto parallel_wall_times =
    parallel_timer().get_summary_data(dealii::TimerOutput::OutputData::total_wall_time);

  for (const auto &[section, wall_time] : serial_wall_times)
    {
      Assert(serial_n_calls.find(section) != serial_n_calls.end(),
             dealii::ExcInternalError());

      const auto n_calls = serial_n_calls.at(section);

      conditionalOStreams::pout_base()
        << "Serial: " << section << " - " << wall_time << "s (" << wall_time / n_calls
        << "s/call, " << n_calls << " calls)"
        << "\n"
        << std::flush;
    }
}

PRISMS_PF_END_NAMESPACE