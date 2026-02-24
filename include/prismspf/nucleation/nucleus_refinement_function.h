// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>
#include <deal.II/grid/tria_accessor.h>

#include <prismspf/core/cell_marker_base.h>
#include <prismspf/core/types.h>

#include <prismspf/nucleation/nucleus.h>

#include <prismspf/user_inputs/nucleation_parameters.h>
#include <prismspf/user_inputs/temporal_discretization.h>

#include <prismspf/config.h>

#include <iostream>
#include <memory>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief The class handles the stochastic nucleation
 * in PRISMS-PF.
 */
template <unsigned int dim>
class NucleusRefinementFunction : public CellMarkerBase<dim>
{
public:
  using CellIterator =
    dealii::CellAccessor<dim>; // dealii::TriaActiveIterator<dealii::CellAccessor<dim>>;
  using prisms::CellMarkerBase<dim>::flag;

  NucleusRefinementFunction(const NucleationParameters      &_nuc_params,
                            const std::vector<Nucleus<dim>> &_nuclei_list)
    : nuc_params(&_nuc_params)
    , nuclei_list(&_nuclei_list)
  {}

  bool
  flag(const CellIterator                     &cell,
       [[maybe_unused]] const SimulationTimer &time_info) const override
  {
    for (const Nucleus<dim> &nucleus : *nuclei_list)
      {
        if (nuc_params->check_active(nucleus, time_info))
          {
            static dealii::Point<dim> unit_corner = []()
            {
              dealii::Point<dim> p;
              for (unsigned int d = 0; d < dim; ++d)
                {
                  p[d] = 1.0;
                }
              return p;
            }();
            dealii::BoundingBox<dim> nucleus_bounding_box(
              std::make_pair<dealii::Point<dim>, dealii::Point<dim>>(
                dealii::Point<dim>(nucleus.location -
                                   (unit_corner * nuc_params->get_refinement_radius())),
                dealii::Point<dim>(nucleus.location +
                                   (unit_corner * nuc_params->get_refinement_radius()))));
            if (cell.bounding_box().has_overlap_with(nucleus_bounding_box))
              {
                return true;
              }
          }
      }
    return false;
  }

  const NucleationParameters      *nuc_params;
  const std::vector<Nucleus<dim>> *nuclei_list;
};

PRISMS_PF_END_NAMESPACE
