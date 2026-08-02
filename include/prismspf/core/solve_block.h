// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/dependencies.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/solve_parameters.h>

#include <prismspf/utilities/assert.h>
#include <prismspf/utilities/logger.h>

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
      ASSERT(
        solve_type == SolveType::Constant || solve_type == SolveType::Explicit ||
          solve_type == SolveType::Linear || solve_type == SolveType::Newton,
        "A valid solve type must be selected (Constant | Explicit | Linear | Newton)",
        solve_type);
      ASSERT(!field_indices.empty(), "Solve blocks must manage at least 1 field");
      if (solve_type == SolveType::Newton)
        {
          for (unsigned int field_index : field_indices)
            {
              const auto &dep_it_rhs = dependencies_rhs.find(field_index);
              ASSERT(dep_it_rhs != dependencies_rhs.end(),
                     "Every field in a newton solve should appear in the residual (RHS) "
                     "expression");
              ASSERT(dep_it_rhs->second.flag != EvalFlags::nothing,
                     "Every field in a newton solve should appear in the residual (RHS) "
                     "expression");
              const auto &dep_it_lhs = dependencies_lhs.find(field_index);
              ASSERT(dep_it_lhs != dependencies_lhs.end(),
                     "Every field in a newton solve should appear as a Delta term in the "
                     "residual Jacobian (LHS) expression");
              ASSERT(dep_it_lhs->second.src_flag != EvalFlags::nothing,
                     "Every field in a newton solve should appear as a Delta term in the "
                     "residual Jacobian (LHS) expression");
            }
        }
      else
        {
          for (unsigned int field_index : field_indices)
            {
              const auto &dep_it_rhs = dependencies_rhs.find(field_index);
              if (dep_it_rhs != dependencies_rhs.end())
                {
                  ASSERT(dep_it_rhs->second.flag == EvalFlags::nothing,
                         "The current value of a field should never appear in the RHS of "
                         "a solve that is not type Newton");
                }
            }
        }
      if (solve_type == SolveType::Linear)
        {
          for (unsigned int field_index : field_indices)
            {
              const auto &dep_it_lhs = dependencies_lhs.find(field_index);
              ASSERT(dep_it_lhs != dependencies_lhs.end(),
                     "Every field in a linear solve should appear in the (LHS) "
                     "expression. Be sure to use the src_flag");
              ASSERT(dep_it_lhs->second.src_flag != EvalFlags::nothing,
                     "Every field in a linear solve should appear in the (LHS) "
                     "expression. Be sure to use the src_flag");
            }
        }
      else if (solve_type == SolveType::Explicit)
        {
          ASSERT(
            dependencies_lhs.empty(),
            "Explicit solves do not have an LHS, and should have no LHS dependencies");
        }
      else if (solve_type == SolveType::Constant)
        {
          ASSERT(dependencies_rhs.empty() && dependencies_lhs.empty(),
                 "Constant \"solves\" do not have an RHS or LHS, and should have no "
                 "dependencies");
        }
      for (const auto &[field_index, dependency] : dependencies_rhs)
        {
          ASSERT(dependency.src_flag == EvalFlags::nothing,
                 "Trial/Change terms should not appear in RHS expressions");
        }
    }
  catch (...)
    {
      Logger::instance() << LogFormatter::error(
                              "Error found during validation of solve block with id " +
                              std::to_string(id))
                         << std::endl;
      throw;
    }
}

inline const std::vector<SolveBlock> &
validate_solve_blocks(const std::vector<SolveBlock>      &solve_blocks,
                      const std::vector<FieldAttributes> &field_attributes)
{
  // Check that field names are unique
  {
    std::set<std::string> field_names;
    for (const auto &field_attribute : field_attributes)
      {
        ASSERT(field_names.find(field_attribute.name) == field_names.end(),
               "Each field must have a unique name");
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
        ASSERT(ids.find(solve_block.id) == ids.end(),
               "Each solve block must have a unique id");
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
            ASSERT(erased == 1,
                   "Each field should be managed by exactly one solve block. The field "
                   "with index " +
                     std::to_string(field_index) +
                     " is anomalously assigned to solve block " +
                     std::to_string(solve_block.id));
          }
      }
    std::string remaining_fields;
    for (unsigned int field_index : field_indices)
      {
        remaining_fields.append(std::to_string(field_index) + ": " +
                                field_attributes[field_index].name + "\n");
      }
    ASSERT(field_indices.empty(),
           "Every field should be managed by exactly one solve block. The following "
           "fields are not managed by any solve block:\n" +
             remaining_fields);
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
            ASSERT(valid_field_indices.find(field_index) != valid_field_indices.end(),
                   "The RHS dependencies of solve block with id " +
                     std::to_string(solve_block.id) +
                     " refer to a field index that is out of range");
          }
        for (const auto &[field_index, dependency] : solve_block.dependencies_lhs)
          {
            ASSERT(valid_field_indices.find(field_index) != valid_field_indices.end(),
                   "The LHS dependencies of solve block with id " +
                     std::to_string(solve_block.id) +
                     " refer to a field index that is out of range");
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
          if (dependency.flag != EvalFlags::nothing &&
              solve_block.solve_type != SolveType::Newton)
            {
              ASSERT(solved_field_indices.find(field_index) != solved_field_indices.end(),
                     "Solve block with id " + std::to_string(solve_block.id) +
                       " has an rhs dependency on the current value of field \"" +
                       field_attributes[field_index].name + "\" with index " +
                       std::to_string(field_index) +
                       " which is not solved in a previous solve block. This is not "
                       "allowed");
            }
        }
      for (const auto &[field_index, dependency] : solve_block.dependencies_lhs)
        {
          if (dependency.flag != EvalFlags::nothing &&
              solve_block.solve_type != SolveType::Newton)
            {
              ASSERT(solved_field_indices.find(field_index) != solved_field_indices.end(),
                     "Solve block with id " + std::to_string(solve_block.id) +
                       " has a lhs dependency on the current value of field \"" +
                       field_attributes[field_index].name + "\" with index " +
                       std::to_string(field_index) +
                       " which is not solved in a previous solve block. This is not "
                       "allowed");
            }
        }
      for (unsigned int field_index : solve_block.field_indices)
        {
          solved_field_indices.insert(field_index);
        }
    }
  return solve_blocks;
}

PRISMS_PF_END_NAMESPACE
