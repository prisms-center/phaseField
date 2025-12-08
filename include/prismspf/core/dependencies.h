// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <array>
#include <set>
#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

namespace Dependencies
{
  static const std::array<EvalFlags, Numbers::max_saved_increments> nothing_arr = []()
  {
    std::array<EvalFlags, Numbers::max_saved_increments> arr {};
    arr.fill(EvalFlags::nothing);
    return arr;
  }();

  // NOLINTBEGIN(misc-non-private-member-variables-in-classes, hicpp-explicit-conversions)
  struct Dependency
  {
    using EvalFlags       = dealii::EvaluationFlags::EvaluationFlags;
    EvalFlags flag        = EvalFlags::nothing;
    EvalFlags change_flag = EvalFlags::nothing;

    std::array<EvalFlags, Numbers::max_saved_increments> old_flags {};

    // Intentional implicit conversions
    /**
     * @brief Construct a Dependency with given flags.
     */
    Dependency(
      EvalFlags _flag        = EvalFlags::nothing,
      EvalFlags _change_flag = EvalFlags::nothing,
      std::array<EvalFlags, Numbers::max_saved_increments> _old_flags = nothing_arr)
      : flag(_flag)
      , change_flag(_change_flag)
      , old_flags(_old_flags)
    {}

    Dependency
    operator|(const Dependency &other) const
    {
      Dependency result(static_cast<EvalFlags>(static_cast<unsigned int>(flag) |
                                               static_cast<unsigned int>(other.flag)),
                        static_cast<EvalFlags>(
                          static_cast<unsigned int>(change_flag) |
                          static_cast<unsigned int>(other.change_flag)));
      for (unsigned int i = 0; i < Numbers::max_saved_increments; ++i)
        {
          result.old_flags.at(i) =
            static_cast<EvalFlags>(static_cast<unsigned int>(old_flags.at(i)) |
                                   static_cast<unsigned int>(other.old_flags.at(i)));
        }
      return result;
    }

    Dependency
    operator&(const Dependency &other) const
    {
      Dependency result(static_cast<EvalFlags>(static_cast<unsigned int>(flag) &
                                               static_cast<unsigned int>(other.flag)),
                        static_cast<EvalFlags>(
                          static_cast<unsigned int>(change_flag) &
                          static_cast<unsigned int>(other.change_flag)));
      for (unsigned int i = 0; i < Numbers::max_saved_increments; ++i)
        {
          result.old_flags.at(i) =
            static_cast<EvalFlags>(static_cast<unsigned int>(old_flags.at(i)) &
                                   static_cast<unsigned int>(other.old_flags.at(i)));
        }
      return result;
    }

    Dependency &
    operator|=(const Dependency &other)
    {
      flag        = static_cast<EvalFlags>(static_cast<unsigned int>(flag) |
                                    static_cast<unsigned int>(other.flag));
      change_flag = static_cast<EvalFlags>(static_cast<unsigned int>(change_flag) |
                                           static_cast<unsigned int>(other.change_flag));
      for (unsigned int i = 0; i < Numbers::max_saved_increments; ++i)
        {
          old_flags.at(i) =
            static_cast<EvalFlags>(static_cast<unsigned int>(old_flags.at(i)) |
                                   static_cast<unsigned int>(other.old_flags.at(i)));
        }

      return *this;
    }

    Dependency &
    operator&=(const Dependency &other)
    {
      flag        = static_cast<EvalFlags>(static_cast<unsigned int>(flag) &
                                    static_cast<unsigned int>(other.flag));
      change_flag = static_cast<EvalFlags>(static_cast<unsigned int>(change_flag) &
                                           static_cast<unsigned int>(other.change_flag));
      for (unsigned int i = 0; i < Numbers::max_saved_increments; ++i)
        {
          old_flags.at(i) =
            static_cast<EvalFlags>(static_cast<unsigned int>(old_flags.at(i)) &
                                   static_cast<unsigned int>(other.old_flags.at(i)));
        }

      return *this;
    }
  };

  // NOLINTEND(misc-non-private-member-variables-in-classes, hicpp-explicit-conversions)

  static const Dependency value {EvalFlags::values};
  static const Dependency gradient {EvalFlags::gradients};
  static const Dependency value_and_gradient {EvalFlags::values | EvalFlags::gradients};
} // namespace Dependencies

using Dependencies::Dependency;

using DependencySet = std::map<Types::Index, Dependency>;

DependencySet
make_dependency_set(const std::vector<FieldAttributes> &field_attributes,
                    const std::set<std::string>        &dependency_strings)
{
  static const std::set<std::pair<std::pair<std::string, std::string>, EvalFlags>>
    reg_delimiters = {
      {{"", ""},           EvalFlags::values   },
      {{"grad(", ")"},     EvalFlags::gradients},
      {{"hess(", ")"},     EvalFlags::hessians },
      {{"hessdiag(", ")"}, EvalFlags::hessians },
      {{"lap(", ")"},      EvalFlags::hessians },
      {{"div(", ")"},      EvalFlags::gradients},
      {{"symgrad(", ")"},  EvalFlags::gradients},
      {{"curl(", ")"},     EvalFlags::gradients},
  };

  DependencySet result;
  for (unsigned int i = 0; i < field_attributes.size(); ++i)
    {
      const auto &attr = field_attributes[i];
      for (const auto &[delimiter, flag] : reg_delimiters)
        {
          std::string potential_match = delimiter.first + attr.name + delimiter.second;
          if (dependency_strings.contains(potential_match))
            {
              result[i].flag |= flag;
            }
          potential_match =
            delimiter.first + "change(" + attr.name + ")" + delimiter.second;
          if (dependency_strings.contains(potential_match))
            {
              result[i].change_flag |= flag;
            }
          for (unsigned int old_index = 0; old_index < Numbers::max_saved_increments;
               ++old_index)
            {
              potential_match = delimiter.first + "old_" + std::to_string(old_index + 1) +
                                "(" + attr.name + ")" + delimiter.second;
              if (dependency_strings.contains(potential_match))
                {
                  result[i].old_flags.at(old_index) |= flag;
                }
            }
        }
    }
  return result;
}

DependencySet
make_dependency_set(const std::vector<FieldAttributes> &field_attributes,
                    const std::string                  &dependency_string)
{
  auto dependency_strings = dealii::Utilities::split_string_list(dependency_string);
  return make_dependency_set(field_attributes,
                             std::set<std::string>(dependency_strings.begin(),
                                                   dependency_strings.end()));
}

PRISMS_PF_END_NAMESPACE