// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/timer.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/timer.h>

#include <prismspf/config.h>

#include <iostream>
#include <mpi.h>
#include <ostream>

#ifdef PRISMS_PF_WITH_CALIPER
#  include <caliper/cali.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

void
Timer::start_section(const char *name)
{
#ifdef PRISMS_PF_WITH_CALIPER
  CALI_MARK_BEGIN(name);
#else
  serial_timer().enter_subsection(name);
#endif
}

void
Timer::end_section(const char *name)
{
#ifdef PRISMS_PF_WITH_CALIPER
  CALI_MARK_END(name);
#else
  serial_timer().leave_subsection();
#endif
}

dealii::TimerOutput &
Timer::serial_timer()
{
  static dealii::TimerOutput instance(ConditionalOStreams::pout_base(),
                                      dealii::TimerOutput::never,
                                      dealii::TimerOutput::wall_times);

  return instance;
}

dealii::TimerOutput &
Timer::parallel_timer()
{
  static dealii::TimerOutput instance(MPI_COMM_WORLD,
                                      std::cout,
                                      dealii::TimerOutput::never,
                                      dealii::TimerOutput::wall_times);

  return instance;
}

void
Timer::print_summary()
{
  // Caliper already prints the summary
#ifdef PRISMS_PF_WITH_CALIPER
  return;
#endif

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
      Assert(serial_n_calls.contains(section), dealii::ExcInternalError());

      const auto n_calls = serial_n_calls.at(section);

      ConditionalOStreams::pout_base()
        << "Serial: " << section << " - " << wall_time << "s (" << wall_time / n_calls
        << "s/call, " << n_calls << " calls)" << "\n"
        << std::flush;
    }
}

PRISMS_PF_END_NAMESPACE