// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

#include <algorithm>

PRISMS_PF_BEGIN_NAMESPACE

struct DependencyExtent
{
  Dependency   dependency;
  unsigned int mg_depth = 0;

  [[nodiscard]] unsigned int
  age() const
  {
    for (unsigned int i = Numbers::max_saved_increments - 1; i >= 0; --i)
      {
        if (dependency.old_flags.at(i) != dealii::EvaluationFlags::nothing)
          {
            return i;
          }
      }
    return 0;
  }

  static std::vector<DependencyExtent>
  calculate_extents(const std::vector<FieldAttributes> &attributes_list,
                    const std::set<SolveGroup>         &solve_groups)
  {
    std::vector<DependencyExtent> extents(attributes_list.size());
    for (const auto &solve_group : solve_groups)
      {
        for (const auto &field_index : solve_group.field_indices)
          {
            DependencyExtent &extent = extents[field_index];
            for (const auto &[index, dep] : solve_group.dependencies_rhs)
              {
                extent.dependency |= dep;
                extent.mg_depth = std::max(extent.mg_depth, solve_group.mg_depth);
              }
          }
      }
    return extents;
  };
};

PRISMS_PF_END_NAMESPACE