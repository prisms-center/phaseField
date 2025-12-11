// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/config.h>
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
Timer::end_section([[maybe_unused]] const char *name)
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
                                      dealii::TimerOutput::summary,
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

  // TODO: Revist the deal.II TimerOutput class and reorganize sections so they are tiered
  // accordingly
}

PRISMS_PF_END_NAMESPACE
