// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

#include <map>
#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE
using EvalFlags = dealii::EvaluationFlags::EvaluationFlags;
using Rank      = FieldInfo::TensorRank;

// NOLINTBEGIN(misc-non-private-member-variables-in-classes)
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
    Rank                      _field_type                  = Rank::Undefined,
    EvalFlags                 _eval_flags_rhs              = EvalFlags::nothing,
    EvalFlags                 _eval_flags_lhs              = EvalFlags::nothing,
    bool                      _is_nucleation_rate_variable = false,
    std::vector<Types::Index> _nucleating_field_indices    = std::vector<Types::Index>())
    : name(std::move(_name))
    , field_type(_field_type)
    , eval_flags_rhs(_eval_flags_rhs)
    , eval_flags_lhs(_eval_flags_lhs)
    , is_nucleation_rate_variable(_is_nucleation_rate_variable)
    , nucleating_field_indices(std::move(_nucleating_field_indices))
  {}

  /**
   * @brief Field name.
   * @remark User-set
   */
  std::string name;

  /**
   * @brief Field type (Scalar/Vector).
   * @remark User-set
   */
  Rank field_type = Rank::Undefined;

  /**
   * @brief Evaluation flags for the types of residual the user is expected to submit to
   * on the RHS.
   * @remark Internally determined
   */
  EvalFlags eval_flags_rhs = EvalFlags::nothing;

  /**
   * @brief Evaluation flags for the types of residual the user is expected to submit to
   * on the LHS. This is empty for Explicit fields.
   * @remark Internally determined
   */
  EvalFlags eval_flags_lhs = EvalFlags::nothing;

  /**
   * @brief Is a nucleation rate
   * @remark User-set
   */
  bool is_nucleation_rate_variable = false;

  /**
   * @brief If this is a nucleation rate, the indices of the nucleating fields
   * @remark User-set
   */
  std::vector<Types::Index> nucleating_field_indices;
};

// NOLINTEND(misc-non-private-member-variables-in-classes)

// TODO: Consider making these static members instead of prismspf::

/**
 * @brief Make a map that maps field names to field indices.
 */
std::map<std::string, Types::Index>
field_index_map(const std::vector<FieldAttributes> &fields)
{
  std::map<std::string, Types::Index> map;
  for (unsigned int i = 0; i < fields.size(); ++i)
    {
      AssertThrow(map.find(fields[i].name) == map.end(),
                  "The names of the fields are not unique. This is not allowed.");
      map[fields[i].name] = i;
    }
  return map;
}

/**
 * @brief Make a map that maps field names to field attributes.
 */
std::map<std::string, FieldAttributes>
field_map(const std::vector<FieldAttributes> &fields)
{
  std::map<std::string, FieldAttributes> map;
  for (const FieldAttributes &field : fields)
    {
      AssertThrow(map.find(field.name) == map.end(),
                  "The names of the fields are not unique. This is not allowed.");
      map[field.name] = field;
    }
  return map;
}

// TODO: Submit a PR/issue to dealii to make operator| const.
// constexpr EvalFlags values_and_gradients = EvalFlags::values | EvalFlags::gradients;

PRISMS_PF_END_NAMESPACE
