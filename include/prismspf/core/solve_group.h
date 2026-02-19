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

  explicit SolveGroup(int                    _id                = -1,
                      PDEType                _pde_type          = PDEType::Explicit,
                      bool                   _is_time_dependent = true,
                      std::set<Types::Index> _field_indices     = {},
                      DependencySet          _dependencies_rhs  = {},
                      DependencySet          _dependencies_lhs  = {})
    : id(_id)
    , pde_type(_pde_type)
    , is_time_dependent(_is_time_dependent)
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
  PDEType pde_type;

  /**
   * @brief Whether the solve is time dependent. This is used to determine whether to
   * initialize the solution vector with the initial conditions or just solve.
   */
  bool is_time_dependent;

  /**
   * @brief Indices of the fields to be solved in this group.
   */
  std::set<Types::Index> field_indices;

  /**
   * @brief Dependencies for the rhs equation(s)
   */
  DependencySet dependencies_rhs;
  /**
   * @brief Dependencies for the lhs equation(s)
   */
  DependencySet dependencies_lhs;

private:
  /**
   * @brief A std::vector used to manage dynamic memory for holding an auxiliary solve
   * group
   * @remark Privately managed
   */
  std::vector<SolveGroup> aux_solve_container;

public:
  void
  set_auxiliary_solve(const SolveGroup &aux_solve = SolveGroup())
  {
    aux_solve_container.clear();
    aux_solve_container.push_back(aux_solve);
  }

  void
  remove_auxiliary_solve() // unset_auxiliary_solve()
  {
    aux_solve_container.clear();
  }

  [[nodiscard]] bool
  has_auxiliary_solve() const
  {
    return !aux_solve_container.empty();
  }

  SolveGroup &
  auxiliary_solve() // get_auxiliary_solve()
  {
    AssertThrow(has_auxiliary_solve(),
                dealii::ExcMessage("Cannot access auxiliary solve when none is set.\n "
                                   "Attempted access from solve id: " +
                                   std::to_string(id)));
    return aux_solve_container[0];
  }

  bool
  operator<(const SolveGroup &other) const
  {
    // TODO: Ensure that the pde_type enum is defined in the proper order.
    // Swap Explicit and ImplicitTimeDependent?
    // Constant | ExplicitTimeDependent | Explicit | ImplicitTimeDependent | Implicit |
    // Postprocess
    if (pde_type < other.pde_type)
      {
        return true;
      }
    return id < other.id;
  }
};

// NOLINTEND(misc-non-private-member-variables-in-classes, hicpp-explicit-conversions)

PRISMS_PF_END_NAMESPACE
