// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/dependencies.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE
inline unsigned int
oldest(Dependencies::Dependency dependencies)
{
  for (unsigned int age = Numbers::max_saved_increments - 1; age >= 0; age--)
    {
      if (dependencies.old_flags.at(age))
        {
          return age;
        }
    }
  return 0;
}

/**
 * @brief Information about what fields need to be held onto. This will likely get
 * refactored to be an oldest age at each level instead of something global like this.
 */
struct DependencyExtents
{
  unsigned int oldest_age   = 0;
  unsigned int max_mg_level = 0;

  // TODO: update this to account multigrid levels. Will need more
  // than std::list<DependencySet>
  DependencyExtents(const std::set<unsigned int>   &field_indices,
                    const std::list<DependencySet> &all_dependeny_sets)
  {
    for (const DependencySet &dependency_set : all_dependeny_sets)
      {
        for (unsigned int field_index : field_indices)
          {
            const auto &dep_pair_it = dependency_set.find(field_index);
            if (dep_pair_it != dependency_set.end())
              {
                oldest_age = std::max(oldest_age, oldest(dep_pair_it->second));
              }
          }
      }
  }
};

PRISMS_PF_END_NAMESPACE