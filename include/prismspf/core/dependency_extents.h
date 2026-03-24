// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/dependencies.h>
#include <prismspf/core/solve_block.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <vector>


PRISMS_PF_BEGIN_NAMESPACE
inline unsigned int
oldest(Dependency dependencies)
{
  for (long int age_index = long(dependencies.old_flags.size()) - 1; age_index >= 0;
       age_index--)
    {
      if (dependencies.old_flags[age_index])
        {
          return age_index + 1;
        }
    }
  return 0;
}

inline int
oldest2(Dependency dependencies)
{
  for (long int age_index = long(dependencies.old_flags.size()) - 1; age_index >= 0;
       age_index--)
    {
      if (dependencies.old_flags[age_index])
        {
          return age_index + 1;
        }
    }
  if (dependencies.flag)
    {
      return 0;
    }
  return -1;
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
  // than std::list<DependencyMap>
  DependencyExtents(const std::set<unsigned int>   &field_indices,
                    const std::list<DependencyMap> &all_dependeny_sets)
  {
    for (const DependencyMap &dependency_set : all_dependeny_sets)
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

struct NewDependencyExtents
{
  std::vector<unsigned int> max_age_per_level;

  NewDependencyExtents(const std::set<unsigned int> &field_indices,
                       const std::list<SolveBlock>  &solve_blocks)
  {
    for (const SolveBlock &solve_block : solve_blocks)
      {
        const unsigned int max_level = 0; // solve_block.max_level; //TODO
        for (unsigned int field_index : field_indices)
          {
            {
              const auto &dep_pair_it = solve_block.dependencies_rhs.find(field_index);
              if (dep_pair_it != solve_block.dependencies_rhs.end())
                {
                  unsigned int age = oldest2(dep_pair_it->second);
                  if (age == -1)
                    {
                      continue;
                    }
                  if (max_level >= max_age_per_level.size())
                    {
                      max_age_per_level.resize(max_level + 1, 0);
                    }
                  for (unsigned int level = 0; level <= max_level; ++level)
                    {
                      max_age_per_level[level] = std::max(max_age_per_level[level], age);
                    }
                }
            }
            {
              const auto &dep_pair_it = solve_block.dependencies_lhs.find(field_index);
              if (dep_pair_it != solve_block.dependencies_lhs.end())
                {
                  unsigned int age = oldest2(dep_pair_it->second);
                  if (age == -1)
                    {
                      continue;
                    }
                  if (max_level >= max_age_per_level.size())
                    {
                      max_age_per_level.resize(max_level + 1, 0);
                    }
                  for (unsigned int level = 0; level <= max_level; ++level)
                    {
                      max_age_per_level[level] = std::max(max_age_per_level[level], age);
                    }
                }
            }
          }
      }
  }
};

PRISMS_PF_END_NAMESPACE
