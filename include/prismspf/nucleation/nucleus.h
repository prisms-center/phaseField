// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/point.h>

#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <ostream>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class contains mutable utilities for phase field problems.
 */
template <unsigned int dim>
struct Nucleus
{
public:
  /**
   * @brief Constructor.
   */
  Nucleus() = default;

  Nucleus(const unsigned int       &_field_index,
          const dealii::Point<dim> &_location,
          const double             &_seed_time,
          const unsigned int       &_seed_increment)
    : field_index(_field_index)
    , location(_location)
    , seed_time(_seed_time)
    , seed_increment(_seed_increment)
  {}

  unsigned int       field_index = 0;
  dealii::Point<dim> location;
  double             seed_time      = 0.0;
  unsigned int       seed_increment = 0;
};

PRISMS_PF_END_NAMESPACE