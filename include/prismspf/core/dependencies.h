// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/conditional_ostreams.h>
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
  // NOLINTBEGIN(misc-non-private-member-variables-in-classes, hicpp-explicit-conversions)
  struct Dependency
  {
    using EvalFlags       = dealii::EvaluationFlags::EvaluationFlags;
    EvalFlags flag        = EvalFlags::nothing;
    EvalFlags change_flag = EvalFlags::nothing;

    std::vector<EvalFlags> old_flags;

    // Intentional implicit conversions
    /**
     * @brief Construct a Dependency with given flags.
     */
    Dependency(EvalFlags                     _flag        = EvalFlags::nothing,
               EvalFlags                     _change_flag = EvalFlags::nothing,
               const std::vector<EvalFlags> &_old_flags   = {})
      : flag(_flag)
      , change_flag(_change_flag)
      , old_flags(_old_flags)
    {}

    Dependency
    operator|(const Dependency &other) const
    {
      Dependency result(flag | other.flag, change_flag | other.change_flag);
      for (unsigned int i = 0; i < std::max(old_flags.size(), other.old_flags.size());
           ++i)
        {
          result.old_flags.push_back(
            (i < old_flags.size() ? old_flags.at(i) : EvalFlags::nothing) |
            (i < other.old_flags.size() ? other.old_flags.at(i) : EvalFlags::nothing));
        }
      return result;
    }

    Dependency
    operator&(const Dependency &other) const
    {
      Dependency result(flag & other.flag, change_flag & other.change_flag);
      for (unsigned int i = 0; i < std::min(old_flags.size(), other.old_flags.size());
           ++i)
        {
          result.old_flags.push_back(old_flags.at(i) & (other.old_flags.at(i)));
        }
      return result;
    }

    Dependency &
    operator|=(const Dependency &other)
    {
      flag        = flag | other.flag;
      change_flag = change_flag | other.change_flag;
      if (other.old_flags.size() > old_flags.size())
        {
          old_flags.resize(other.old_flags.size());
        }
      for (unsigned int i = 0; i < old_flags.size(); ++i)
        {
          old_flags.at(i) = old_flags.at(i) | other.old_flags.at(i);
        }
      return *this;
    }

    Dependency &
    operator&=(const Dependency &other)
    {
      flag        = flag & other.flag;
      change_flag = change_flag & other.change_flag;
      old_flags.resize(std::min(old_flags.size(), other.old_flags.size()));
      for (unsigned int i = 0; i < old_flags.size(); ++i)
        {
          old_flags.at(i) = old_flags.at(i) & other.old_flags.at(i);
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

constexpr unsigned int default_max_checked_age = 6;

template <unsigned int max_checked_age = default_max_checked_age>
inline DependencySet
make_dependency_set(const std::vector<FieldAttributes> &field_attributes,
                    std::set<std::string>               dependency_strings)
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
          auto        iter            = dependency_strings.find(potential_match);
          if (iter != dependency_strings.end())
            {
              result[i].flag |= flag;
              dependency_strings.erase(iter);
            }
          potential_match =
            delimiter.first + "change(" + attr.name + ")" + delimiter.second;
          iter = dependency_strings.find(potential_match);
          if (iter != dependency_strings.end())
            {
              result[i].change_flag |= flag;
              dependency_strings.erase(iter);
            }
          for (unsigned int old_index = 0; old_index < max_checked_age; ++old_index)
            {
              potential_match = delimiter.first + "old_" + std::to_string(old_index + 1) +
                                "(" + attr.name + ")" + delimiter.second;
              iter = dependency_strings.find(potential_match);
              if (iter != dependency_strings.end())
                {
                  if (result[i].old_flags.size() < old_index + 1)
                    {
                      result[i].old_flags.resize(old_index + 1, EvalFlags::nothing);
                    }
                  result[i].old_flags.at(old_index) |= flag;
                  dependency_strings.erase(iter);
                }
            }
        }
    }
  if (!dependency_strings.empty())
    {
      std::string error_message = "The following dependencies were not recognized: ";
      for (const auto &dep : dependency_strings)
        {
          error_message += dep + ", ";
        }
      error_message.pop_back();
      error_message.pop_back();
#ifndef DEBUG
      ConditionalOStreams::pout_base() << "Warning: " << error_message << std::endl;
#endif
      Assert(false, dealii::ExcMessage(error_message));
    }
  return result;
}

inline DependencySet
make_dependency_set(const std::vector<FieldAttributes> &field_attributes,
                    const std::string                  &dependency_string)
{
  auto dependency_strings = dealii::Utilities::split_string_list(dependency_string);
  return make_dependency_set(field_attributes,
                             std::set<std::string>(dependency_strings.begin(),
                                                   dependency_strings.end()));
}

PRISMS_PF_END_NAMESPACE