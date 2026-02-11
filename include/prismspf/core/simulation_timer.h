// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/types.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

class SimulationTimer
{
public:
  SimulationTimer() = default;

  explicit SimulationTimer(double step_size)
    : time_step_size(step_size)
  {}

  [[nodiscard]] unsigned int
  get_increment() const
  {
    return current_increment;
  }

  [[nodiscard]] double
  get_time() const
  {
    return current_time;
  }

  [[nodiscard]] double
  get_timestep() const
  {
    return time_step_size;
  }

  void
  increment(double step_size)
  {
    current_increment++;
    current_time += step_size;
  }

  void
  increment()
  {
    increment(time_step_size);
  }

  void
  set_timestep(double step_size)
  {
    time_step_size = step_size;
  }

  void
  set_time(double time)
  {
    current_time = time;
  }

  void
  set_increment(unsigned int increment)
  {
    current_increment = increment;
  }

  void
  reset()
  {
    current_increment = 0;
    current_time      = 0.0;
  }

private:
  unsigned int current_increment = 0;
  double       current_time      = 0.0;
  double       time_step_size    = 0.0;
};

PRISMS_PF_END_NAMESPACE