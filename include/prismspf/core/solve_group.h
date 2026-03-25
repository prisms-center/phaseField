// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/dependencies.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <set>
#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Enum describing when each group of fields gets solved.
 */
enum SolveTiming
{
  /**
   * @brief Primary fields are initialized explicitly through initial conditions rather
   * than through the solver on increment zero.
   * @example Concentration variable
   */
  Primary,
  Initialized = Primary,
  /**
   * @brief Secondary fields are only evaluated by the pde solver on every increment, not
   * initialized by a separate function.
   * @example Chemical potential
   */
  Secondary,
  Uninitialized = Secondary,
  /**
   * @brief PostProcess fields are only solved on output increments.
   * @example Local free energy density
   */
  PostProcess,
  /**
   * @brief NucleationRate fields are only solved on nucleation attempt and output
   * increments.
   */
  NucleationRate
};

// NOLINTBEGIN(misc-non-private-member-variables-in-classes, hicpp-explicit-conversions)
// readability-simplify-boolean-expr

/**
 * @brief Structure to hold the attributes of a solve-group.
 */
class SolveGroup
{
public:
  using EvalFlags = dealii::EvaluationFlags::EvaluationFlags;
  using FieldType = TensorRank;

  explicit SolveGroup(int                    _id               = -1,
                      SolveType              _solve_type       = Explicit,
                      SolveTiming            _solve_timing     = Primary,
                      std::set<Types::Index> _field_indices    = {},
                      DependencyMap          _dependencies_rhs = {},
                      DependencyMap          _dependencies_lhs = {})
    : id(_id)
    , solve_type(_solve_type)
    , solve_timing(_solve_timing)
    , field_indices(std::move(_field_indices))
    , dependencies_rhs(std::move(_dependencies_rhs))
    , dependencies_lhs(std::move(_dependencies_lhs))
  {}

  /**
   * @brief Unique identifier. Use this in 'if' statements or switch cases in equations
   * lhs and rhs.
   */
  int id;

  /**
   * @brief PDE type (Constant | Explicit | Linear | Newton).
   */
  SolveType solve_type;

  /**
   * @brief This is used to determine whether to
   * initialize the solution vector with the initial conditions or just solve.
   */
  SolveTiming solve_timing;

  /**
   * @brief Indices of the fields to be solved in this group.
   */
  std::set<Types::Index> field_indices;

  /**
   * @brief Dependencies for the rhs equation(s)
   */
  DependencyMap dependencies_rhs;
  /**
   * @brief Dependencies for the lhs equation(s)
   */
  DependencyMap dependencies_lhs;

  /**
   * @brief Solves that occur inside the parent solve. This is only really meant for
   * implicit solves with higher order spatial derivatives.
   */
  std::vector<SolveGroup> aux_solve_container;

  [[nodiscard]] bool
  has_auxiliary_solve() const
  {
    return !aux_solve_container.empty();
  }

  bool
  operator<(const SolveGroup &other) const
  {
    return id < other.id;
  }

  void
  validate() const
  {
    AssertThrow(
      solve_type == SolveType::Constant || solve_type == SolveType::Explicit ||
        solve_type == SolveType::Linear || solve_type == SolveType::Newton,
      dealii::ExcMessage(
        "A valid solve type must be selected (Constant | Explicit | Linear | Newton)\n"));
    AssertThrow(!field_indices.empty(),
                dealii::ExcMessage("This solve group must manage at least 1 field.\n"));
    if (solve_type == SolveType::Newton)
      {
        for (unsigned int field_index : field_indices)
          {
            const auto &dep_it_rhs = dependencies_rhs.find(field_index);
            AssertThrow(dep_it_rhs != dependencies_rhs.end(),
                        dealii::ExcMessage("Every field in a newton solve should appear "
                                           "in the residual (RHS) expression.\n"));
            AssertThrow(dep_it_rhs->second.flag != EvalFlags::nothing,
                        dealii::ExcMessage("Every field in a newton solve should appear "
                                           "in the residual (RHS) expression.\n"));
            const auto &dep_it_lhs = dependencies_lhs.find(field_index);
            AssertThrow(dep_it_lhs != dependencies_lhs.end(),
                        dealii::ExcMessage(
                          "Every field in a newton solve should appear as a Delta term"
                          "in the residual Jacobian (LHS) expression.\n"));
            AssertThrow(dep_it_lhs->second.src_flag != EvalFlags::nothing,
                        dealii::ExcMessage(
                          "Every field in a newton solve should appear as a Delta term"
                          "in the residual Jacobian (LHS) expression.\n"));
          }
      }
    else if (solve_type == SolveType::Linear)
      {
        for (unsigned int field_index : field_indices)
          {
            const auto &dep_it_lhs = dependencies_lhs.find(field_index);
            AssertThrow(dep_it_lhs != dependencies_lhs.end(),
                        dealii::ExcMessage(
                          "Every field in a linear solve should appear "
                          "in the (LHS) expression. Be sure to use the src_flag.\n"));
            AssertThrow(dep_it_lhs->second.src_flag != EvalFlags::nothing,
                        dealii::ExcMessage(
                          "Every field in a linear solve should appear "
                          "in the (LHS) expression. Be sure to use the src_flag.\n"));
          }
      }
    else if (solve_type == SolveType::Explicit)
      {
        AssertThrow(dependencies_lhs.empty(),
                    dealii::ExcMessage("Explicit solves do not have an LHS, "
                                       "and should have no LHS dependencies.\n"));
      }
    else if (solve_type == SolveType::Constant)
      {
        AssertThrow(dependencies_rhs.empty() && dependencies_lhs.empty(),
                    dealii::ExcMessage("Constant \"solves\" do not have an RHS or LHS, "
                                       "and should have no dependencies.\n"));
      }
    for (const auto &[field_index, dependency] : dependencies_rhs)
      {
        AssertThrow(dependency.src_flag == EvalFlags::nothing,
                    dealii::ExcMessage(
                      "Trial/Change terms should not appear in RHS expressions.\n"));
      }
    for (const SolveGroup &aux : aux_solve_container)
      {
        aux.validate();
      }
  }
};

// NOLINTEND(misc-non-private-member-variables-in-classes, hicpp-explicit-conversions)

PRISMS_PF_END_NAMESPACE
