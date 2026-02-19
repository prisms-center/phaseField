// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/timer.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Timer class for PRISMS-PF
 *
 * This class keeps track of nested timing sections so we can print a structured summary
 * at the end of the simulation.
 *
 * Fortunately, `Caliper` handles this nicely (and more!) or us. `deal.II` doesn't do this
 * so we have to keep track of a few additional objects. The logical way to represent this
 * is with a tree node structure.
 */
class Timer
{
public:
  /**
   * @brief Constructor.
   */
  Timer() = default;

  /**
   * @brief Destructor.
   *
   * Calls `print_summary` upon destruction.
   */
  ~Timer();

  Timer(const Timer &) = delete;
  Timer &
  operator=(const Timer &) = delete;
  Timer(Timer &&)          = delete;
  Timer &
  operator=(Timer &&) = delete;

  /**
   * @brief Timer scope guard.
   *
   * This allows you to use the timer like the following
   * @code
   * void f()
   * {
   *  Timer::Scope outer("outer");
   *
   *  {
   *    Timer::Scope inner("inner");
   *    // Work goes here
   *  } // Inner scope ends here
   *
   * // Work goes here
   * } // Outer scope ends here
   * @endcode
   *
   */
  struct Scope
  {
  public:
    explicit Scope(const char *name)
      : name(name)
    {
      Timer::start_section(name);
    }

    ~Scope()
    {
      Timer::end_section(name);
    }

    Scope(const Scope &) = delete;
    Scope &
    operator=(const Scope &) = delete;
    Scope(Scope &&)          = delete;
    Scope &
    operator=(Scope &&) = delete;

  private:
    const char *name;
  };

  /**
   * @brief Start a new timer section.
   */
  static void
  start_section(const char *name);

  /**
   * @brief End the timer section.
   */
  static void
  end_section(const char *name);

  /**
   * @brief deal.II timer for the 0th MPI process
   */
  static dealii::TimerOutput &
  serial_timer();

  /**
   * @brief deal.II timer for parallel MPI process
   */
  static dealii::TimerOutput &
  parallel_timer();

  /**
   * @brief Print a sorted summary of the timed sections.
   */
  static void
  print_summary();
};

PRISMS_PF_END_NAMESPACE
