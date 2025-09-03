// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <prismspf/user_inputs/temporal_discretization.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Base class for cell markers.
 *
 * This class provides an interface for marking cells in the simulation domain for
 * refinement based on position.
 */
template <unsigned int dim>
class CellMarkerBase
{
public:
  using CellIterator =
    dealii::CellAccessor<dim>; // dealii::TriaActiveIterator<dealii::CellAccessor<dim>>;
  /**
   * @brief Destructor.
   */
  virtual ~CellMarkerBase() = default;

  virtual bool
  flag([[maybe_unused]] const CellIterator           &cell,
       [[maybe_unused]] const TemporalDiscretization &time_info) const
  {
    return this->flag(cell.center(), time_info);
  }

  virtual bool
  flag([[maybe_unused]] const dealii::Point<dim>     &point,
       [[maybe_unused]] const TemporalDiscretization &time_info) const
  {
    return false;
  }
};

PRISMS_PF_END_NAMESPACE