// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/dependencies.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/linear_solve_parameters.h>
#include <prismspf/user_inputs/nonlinear_solve_parameters.h>

#include <prismspf/config.h>

#include <set>
#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Enum describing when each block of fields gets solved.
 */
enum SolveTiming
{
  /**
   * @brief Primary fields are initialized explicitly through initial conditions rather
   * than through the solver on increment zero.
   */
  Primary,
  Initialized = Primary,
  /**
   * @brief Secondary fields are only evaluated by the pde solver on every increment, not
   * initialized by a separate function.
   */
  Secondary,
  Uninitialized = Secondary,
  /**
   * @brief PostProcess fields are only solved on output increments.
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
 * @brief Structure to hold the attributes of a solve-block.
 */
class SolveBlock
{
public:
  using EvalFlags = dealii::EvaluationFlags::EvaluationFlags;
  using FieldType = TensorRank;

  explicit SolveBlock(int                    _id               = -1,
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
   * @brief Indices of the fields to be solved in this block.
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
   * @brief Linear solver parameters. Only used for linear and newton solve blocks.
   * @note May be overridden by user input parameters.
   */
  LinearSolverParameters linear_solver_parameters;

  /**
   * @brief Linear solver parameters. Only used for linear and newton solve blocks.
   * @note May be overridden by user input parameters.
   */
  NonlinearSolverParameters nonlinear_solver_parameters;

  bool
  operator<(const SolveBlock &other) const
  {
    return id < other.id;
  }

  void
  validate() const;
};

inline void
SolveBlock::validate() const
{
  try
    {
      AssertThrow(solve_type == SolveType::Constant ||
                    solve_type == SolveType::Explicit ||
                    solve_type == SolveType::Linear || solve_type == SolveType::Newton,
                  dealii::ExcMessage("A valid solve type must be selected (Constant | "
                                     "Explicit | Linear | Newton)\n"));
      AssertThrow(!field_indices.empty(),
                  dealii::ExcMessage("This solve block must manage at least 1 field.\n"));
      if (solve_type == SolveType::Newton)
        {
          for (unsigned int field_index : field_indices)
            {
              const auto &dep_it_rhs = dependencies_rhs.find(field_index);
              AssertThrow(dep_it_rhs != dependencies_rhs.end(),
                          dealii::ExcMessage(
                            "Every field in a newton solve should appear "
                            "in the residual (RHS) expression.\n"));
              AssertThrow(dep_it_rhs->second.flag != EvalFlags::nothing,
                          dealii::ExcMessage(
                            "Every field in a newton solve should appear "
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
      else
        {
          for (unsigned int field_index : field_indices)
            {
              const auto &dep_it_rhs = dependencies_rhs.find(field_index);
              if (dep_it_rhs != dependencies_rhs.end())
                {
                  AssertThrow(dep_it_rhs->second.flag == EvalFlags::nothing,
                              dealii::ExcMessage(
                                "The current value of a field should never appear "
                                "in the RHS of a solve that is not type Newton.\n"));
                }
            }
        }
      if (solve_type == SolveType::Linear)
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
    }
  catch (...)
    {
      ConditionalOStreams::pout_base()
        << "Error found during validation of solve block with id " << id << "\n";
      throw;
    }
}

inline void
validate_solve_blocks(const std::vector<SolveBlock>      &solve_blocks,
                      const std::vector<FieldAttributes> &field_attributes)
{
  // Check that field names are unique
  {
    std::set<std::string> field_names;
    for (const auto &field_attribute : field_attributes)
      {
        AssertThrow(field_names.find(field_attribute.name) == field_names.end(),
                    dealii::ExcMessage("Each field must have a unique name.\n"));
        field_names.insert(field_attribute.name);
      }
  }
  // Validate each solve block individually
  for (const auto &solve_block : solve_blocks)
    {
      solve_block.validate();
    }
  // Check for duplicate solve block ids
  {
    std::set<int> ids;
    for (const auto &solve_block : solve_blocks)
      {
        AssertThrow(ids.find(solve_block.id) == ids.end(),
                    dealii::ExcMessage("Each solve block must have a unique id.\n"));
        ids.insert(solve_block.id);
      }
  }
  // Check that each field is in exactly one solve block
  {
    std::set<unsigned int> field_indices;
    for (unsigned int i = 0; i < field_attributes.size(); ++i)
      {
        field_indices.insert(i);
      }
    for (const auto &solve_block : solve_blocks)
      {
        for (unsigned int field_index : solve_block.field_indices)
          {
            size_t erased = field_indices.erase(field_index);
            AssertThrow(erased == 1,
                        dealii::ExcMessage("Each field should be managed by exactly one "
                                           "solve block. The field with index " +
                                           std::to_string(field_index) +
                                           " is anomalously assigned to solve block " +
                                           std::to_string(solve_block.id) + ".\n"));
          }
      }
    std::string remaining_fields;
    for (unsigned int field_index : field_indices)
      {
        remaining_fields.append(std::to_string(field_index) + ": " +
                                field_attributes[field_index].name + "\n");
      }
    AssertThrow(
      field_indices.empty(),
      dealii::ExcMessage(
        "Every field should be managed by exactly one solve block. The following "
        "fields are not managed by any solve block:\n" +
        remaining_fields));
  }
  // Check that the dependencies of each solve block only refer to fields that are
  // defined in the field attributes.
  {
    std::set<unsigned int> valid_field_indices;
    for (unsigned int i = 0; i < field_attributes.size(); ++i)
      {
        valid_field_indices.insert(i);
      }
    for (const auto &solve_block : solve_blocks)
      {
        for (const auto &[field_index, dependency] : solve_block.dependencies_rhs)
          {
            AssertThrow(
              valid_field_indices.find(field_index) != valid_field_indices.end(),
              dealii::ExcMessage("The RHS dependencies of solve block with id " +
                                 std::to_string(solve_block.id) +
                                 " refer to a field index that is out of range.\n"));
          }
        for (const auto &[field_index, dependency] : solve_block.dependencies_lhs)
          {
            AssertThrow(
              valid_field_indices.find(field_index) != valid_field_indices.end(),
              dealii::ExcMessage("The LHS dependencies of solve block with id " +
                                 std::to_string(solve_block.id) +
                                 " refer to a field index that is out of range.\n"));
          }
      }
  }
  // Check that the order of the solve blocks is consistent with their solve timing.
  // This will be non-exhaustive. Only checking that there are no dependencies for current
  // value that doesn't exist yet.

  std::set<unsigned int> solved_field_indices;
  for (const auto &solve_block : solve_blocks)
    {
      for (const auto &[field_index, dependency] : solve_block.dependencies_rhs)
        {
          if (dependency.flag != EvalFlags::nothing)
            {
              AssertThrow(
                solved_field_indices.find(field_index) != solved_field_indices.end(),
                dealii::ExcMessage(
                  "Solve block with id " + std::to_string(solve_block.id) +
                  " has a dependency on the current value of field with index " +
                  std::to_string(field_index) +
                  " which is not solved in a previous solve block. This is not "
                  "allowed.\n"));
            }
        }
      for (const auto &[field_index, dependency] : solve_block.dependencies_lhs)
        {
          if (dependency.flag != EvalFlags::nothing &&
              solve_block.solve_type != SolveType::Newton)
            {
              AssertThrow(
                solved_field_indices.find(field_index) != solved_field_indices.end(),
                dealii::ExcMessage(
                  "Solve block with id " + std::to_string(solve_block.id) +
                  " has a dependency on the current value of field with index " +
                  std::to_string(field_index) +
                  " which is not solved in a previous solve block. This is not "
                  "allowed.\n"));
            }
        }
      for (unsigned int field_index : solve_block.field_indices)
        {
          solved_field_indices.insert(field_index);
        }
    }
}

using SolveGroup = SolveBlock;

// NOLINTEND(misc-non-private-member-variables-in-classes, hicpp-explicit-conversions)

PRISMS_PF_END_NAMESPACE
