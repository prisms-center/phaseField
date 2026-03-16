// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <map>
#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Structure to hold the attributes of a field. This includes things like
 * the name, rank, and nucleation information..
 */
struct FieldAttributes
{
  /**
   * @brief Constructor
   */
  explicit FieldAttributes(
    std::string               _name                        = "",
    TensorRank                _field_type                  = TensorRank::Scalar,
    bool                      _is_nucleation_rate_variable = false,
    std::vector<Types::Index> _nucleating_field_indices    = std::vector<Types::Index>())
    : name(std::move(_name))
    , field_type(_field_type)
    , is_nucleation_rate_variable(_is_nucleation_rate_variable)
    , nucleating_field_indices(std::move(_nucleating_field_indices))
  {}

  /**
   * @brief Field name.
   */
  std::string name;

  /**
   * @brief Field type (Scalar/Vector).
   */
  TensorRank field_type = TensorRank::Scalar;

  /**
   * @brief Is a nucleation rate
   */
  bool is_nucleation_rate_variable = false;

  /**
   * @brief If this is a nucleation rate, the indices of the nucleating fields
   */
  std::vector<Types::Index> nucleating_field_indices;
};

/**
 * @brief Make a map that maps field names to field indices.
 */
inline std::map<std::string, Types::Index>
field_index_map(const std::vector<FieldAttributes> &fields)
{
  std::map<std::string, Types::Index> map;
  for (unsigned int i = 0; i < fields.size(); ++i)
    {
      AssertThrow(map.find(fields[i].name) == map.end(),
                  dealii::ExcMessage(
                    "The names of the fields are not unique. This is not allowed."));
      map[fields[i].name] = i;
    }
  return map;
}

/**
 * @brief Make a map that maps field names to field attributes.
 */
inline std::map<std::string, FieldAttributes>
field_map(const std::vector<FieldAttributes> &fields)
{
  std::map<std::string, FieldAttributes> map;
  for (const FieldAttributes &field : fields)
    {
      AssertThrow(map.find(field.name) == map.end(),
                  dealii::ExcMessage(
                    "The names of the fields are not unique. This is not allowed."));
      map[field.name] = field;
    }
  return map;
}

// TODO: Submit a PR/issue to dealii to make operator| constexpr.
// constexpr EvalFlags values_and_gradients = EvalFlags::values |
// EvalFlags::gradients;

PRISMS_PF_END_NAMESPACE
