// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/timer.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/timer.h>

#include <prismspf/config.h>

#include <iomanip>
#include <iostream>
#include <map>
#include <mpi.h>
#include <stack>
#include <string>
#include <vector>

#ifdef PRISMS_PF_WITH_CALIPER
#  include <caliper/cali.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

namespace
{
  struct TimerStack
  {
    /**
     * @brief Stack of active section names.
     */
    std::stack<std::string> active;

    /**
     * @brief List of sections, ordered by insertion order, for the summary.
     */
    std::vector<std::string> insertion_order;

    /**
     * @brief Section depth and parent key.
     *
     * Parent key is empty if at the top of hierarchy.
     */
    struct Meta
    {
      unsigned int depth = 0;
      std::string  parent;
    };

    /**
     * @brief Parent and depth of each section.
     */
    std::map<std::string, Meta> meta;

    /**
     * @brief Get the current active section.
     */
    [[nodiscard]] std::string
    current_section() const
    {
      return active.empty() ? "" : active.top();
    }

    /**
     * @brief Add a new section.
     *
     * Returns the key string that should be passed to the deal.II timer class.
     */
    std::string
    push(const std::string &name)
    {
      const std::string parent = current_section();
      const std::string key    = parent.empty() ? name : parent + " > " + name;

      if (!meta.contains(key))
        {
          Meta data;
          data.depth  = static_cast<unsigned int>(active.size());
          data.parent = parent;
          meta[key]   = data;
          insertion_order.push_back(key);
        }

      active.push(key);
      return key;
    }

    /**
     * @brief Remove the current section.
     */
    void
    pop()
    {
      AssertThrow(!active.empty(), dealii::ExcMessage("Timer stack underflow"));
      active.pop();
    }
  };

  TimerStack &
  timer_stack()
  {
    static TimerStack instance;
    return instance;
  }

} // namespace

void
Timer::start_section(const char *name)
{
#ifdef PRISMS_PF_WITH_CALIPER
  CALI_MARK_BEGIN(name);
#else
  const std::string key = timer_stack().push(name);
  serial_timer().enter_subsection(key);
#endif
}

void
Timer::end_section(const char *name)
{
#ifdef PRISMS_PF_WITH_CALIPER
  CALI_MARK_END(name);
#else
  // Check that the name matches the active section
  const std::string &active_section = timer_stack().current_section();
  AssertThrow(
    !active_section.empty() &&
      (active_section == name ||
       active_section.size() > std::string(name).size() &&
         active_section.substr(active_section.size() - std::string(name).size()) == name),
    dealii::ExcMessage(std::string("Timer::end_section mismatch: expected segment '") +
                       name + "' but top of stack is '" + active_section + "'."));
  serial_timer().leave_subsection();
  timer_stack().pop();
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

  const auto &stack = timer_stack();
  const auto  serial_n_calls_data =
    serial_timer().get_summary_data(dealii::TimerOutput::OutputData::n_calls);
  const auto serial_wall_time_data =
    serial_timer().get_summary_data(dealii::TimerOutput::OutputData::total_wall_time);
  const auto serial_cpu_time_data =
    serial_timer().get_summary_data(dealii::TimerOutput::OutputData::total_cpu_time);

  auto lookup = [](const std::map<std::string, double> &data,
                   const std::string                   &key) -> double
    {
      const auto iterator = data.find(key);
      return iterator != data.end() ? iterator->second : 0.0;
    };

  // Column widths
  constexpr int w_label = 44;
  constexpr int w_calls = 8;
  constexpr int w_wall  = 12;
  constexpr int w_cpu   = 12;
  constexpr int w_pct   = 9;
  constexpr int total_w = w_label + w_calls + w_wall + w_cpu + w_pct;

  auto &out = ConditionalOStreams::pout_base();

  out << "\n"
      << std::string(total_w, '=') << "\n"
      << "  PRISMS-PF Timing Summary\n"
      << std::string(total_w, '=') << "\n"
      << std::left << std::setw(w_label) << "Section" << std::right << std::setw(w_calls)
      << "N Calls" << std::setw(w_wall) << "Wall Time (s)" << std::setw(w_cpu)
      << "CPU Time (s)" << std::setw(w_pct) << "% Parent."
      << "\n"
      << std::string(total_w, '-') << "\n";

  for (const auto &key : stack.insertion_order)
    {
      const auto        &meta  = stack.meta.at(key);
      const unsigned int depth = meta.depth;

      const double wall_time = lookup(serial_wall_time_data, key);
      const double cpu_time  = lookup(serial_cpu_time_data, key);
      const auto   n_calls = static_cast<unsigned int>(lookup(serial_n_calls_data, key));

      // Percentage relative to parent (or 100% for roots)
      double pct = 100.0;
      if (!meta.parent.empty())
        {
          const double parent_wall = lookup(serial_wall_time_data, meta.parent);
          if (parent_wall > 0.0)
            {
              pct = wall_time / parent_wall * 100.0;
            }
        }

      // Bare section name (last segment after " > ")
      const std::size_t seperator = key.rfind(" > ");
      const std::string bare =
        seperator == std::string::npos ? key : key.substr(seperator + 3);

      const std::string indent(static_cast<std::size_t>(depth) * 2, ' ');
      const std::string label = indent + (depth > 0 ? "|- " : "") + bare;

      out << std::left << std::setw(w_label) << label << std::right << std::fixed
          << std::setprecision(3) << std::setw(w_calls) << n_calls << std::setw(w_wall)
          << wall_time << std::setw(w_cpu) << cpu_time << std::setw(w_pct - 1) << pct
          << "%"
          << "\n";
    }

  out << std::string(total_w, '=') << "\n\n";
}

PRISMS_PF_END_NAMESPACE
