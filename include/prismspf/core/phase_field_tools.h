// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/types.h>

#include <prismspf/nucleation/nucleus.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class contains mutable utilities for phase field problems.
 */
template <unsigned int dim>
struct PhaseFieldTools
{
  /**
   * @brief Nucleus list.
   */
  std::vector<Nucleus<dim>> nuclei_list;
};

PRISMS_PF_END_NAMESPACE
