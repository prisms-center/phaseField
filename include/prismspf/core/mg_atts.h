// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Utility to keep track of multi-level indices
 */
struct MultiLevelAtts
{
  // Reminder MGLevelObject constructor is *inclusive*, so these are the actual min and
  // max levels we have.
  unsigned int min_level;
  unsigned int max_level;

  [[nodiscard]] unsigned int
  num_levels() const
  {
    return max_level - min_level + 1;
  }

  /**
   * @brief Get the relative level from the global level.
   */
  [[nodiscard]] unsigned int
  get_relative_level(unsigned int level) const
  {
    return max_level - level;
  }

  /**
   * @brief Get the global level from the relative level.
   */
  [[nodiscard]] unsigned int
  get_level(unsigned int relative_level) const
  {
    return max_level - relative_level;
  }
};

PRISMS_PF_END_NAMESPACE