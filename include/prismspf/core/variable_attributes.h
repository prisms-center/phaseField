// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Structure to hold the raw string dependencies of a field.
 */
struct RawDependencies
{
  /**
   * @brief The user-inputted dependencies for the RHS value term.
   * @remark User-set
   */
  std::set<std::string> dependencies_value_rhs;

  /**
   * @brief The user-inputted dependencies for the RHS gradient term.
   * @remark User-set
   */
  std::set<std::string> dependencies_gradient_rhs;

  /**
   * @brief The user-inputted dependencies for the LHS value term.
   * @remark User-set
   */
  std::set<std::string> dependencies_value_lhs;

  /**
   * @brief The user-inputted dependencies for the LHS gradient term.
   * @remark User-set
   */
  std::set<std::string> dependencies_gradient_lhs;

  /**
   * @brief The collection of value and gradient dependencies for the RHS.
   *
   * This is just a superset of `dependencies_value_rhs` and `dependencies_gradient_rhs`.
   * @remark Internally determined
   */
  std::set<std::string> dependencies_rhs;

  /**
   * @brief The collection of value and gradient dependencies for the LHS.
   *
   * This is just a superset of `dependencies_value_lhs` and `dependencies_gradient_lhs`.
   * @remark Internally determined
   */
  std::set<std::string> dependencies_lhs;

  /**
   * @brief The collection of nucleating fields.
   * @remark User-set
   */
  std::set<std::string> nucleating_fields;
};

/**
 * @brief Information about the field.
 */
struct FieldInfo
{
  /**
   * @brief Tensor rank of the field.
   *
   * Currently, only scalar and vectors are supported.
   */
  enum class TensorRank : std::uint8_t
  {
    Undefined = static_cast<std::uint8_t>(-1),
    Scalar    = 0,
    Vector    = 1
  };

  /**
   * @brief Type of PDE that is being solved.
   */
  enum class PDEType : std::uint8_t
  {
    Undefined,
    Constant,
    ExplicitTimeDependent,
    ImplicitTimeDependent,
    TimeIndependent,
    Auxiliary
  };

  /**
   * @brief Internal classification of PDE types.
   *
   * There are several different types of solves that are possible.
   *
   * For PDEType::ExplicitTimeDependent, all fields of that type can be solved
   * concurrently. This also applies for postprocessed fields of that type.
   *
   * For PDEType::Constant, those fields are never updated.
   *
   * For the rest, we have several classifications that can either be solved sequentially
   * of concurrently.
   */
  enum class SolveType : std::uint8_t
  {
    Undefined,
    ConcurrentConstant,
    ConcurrentExplicit,
    ConcurrentAuxiliary,
    SequentialAuxiliary,
    SequentialLinear,
    SequentialSelfNonlinear,
    SequentialCoNonlinear
  };

  /**
   * @brief Internal classification of dependency type flags.
   */
  enum class DependencyTypeFlags : std::uint8_t
  {
    None = 0x00,

    New     = 0x01,
    Current = 0x02,

    Change = 0x04,

    OldOne   = 0x08,
    OldTwo   = 0x10,
    OldThree = 0x20,
    OldFour  = 0x40
  };

  friend DependencyTypeFlags
  operator|(const DependencyTypeFlags flag_1, const DependencyTypeFlags flag_2)
  {
    return static_cast<DependencyTypeFlags>(static_cast<std::uint8_t>(flag_1) |
                                            static_cast<std::uint8_t>(flag_2));
  }

  friend DependencyTypeFlags &
  operator|=(DependencyTypeFlags &flag_1, const DependencyTypeFlags flag_2)
  {
    flag_1 = flag_1 | flag_2;
    return flag_1;
  }

  friend DependencyTypeFlags
  operator&(const DependencyTypeFlags flag_1, const DependencyTypeFlags flag_2)
  {
    return static_cast<DependencyTypeFlags>(static_cast<std::uint8_t>(flag_1) &
                                            static_cast<std::uint8_t>(flag_2));
  }

  friend DependencyTypeFlags &
  operator&=(DependencyTypeFlags &flag_1, const DependencyTypeFlags flag_2)
  {
    flag_1 = flag_1 & flag_2;
    return flag_1;
  }

  bool                is_postprocess        = false;
  bool                is_nucleation_rate    = false;
  TensorRank          tensor_rank           = TensorRank::Undefined;
  PDEType             pde_type              = PDEType::Undefined;
  SolveType           solve_type            = SolveType::Undefined;
  DependencyTypeFlags dependency_type_flags = DependencyTypeFlags::None;
  unsigned int        global_index          = Numbers::invalid_index;
  unsigned int        solve_block_index     = Numbers::invalid_index;
};

inline std::string
to_string(FieldInfo::TensorRank tensor_rank)
{
  switch (tensor_rank)
    {
      case FieldInfo::TensorRank::Undefined:
        return "Undefined";
      case FieldInfo::TensorRank::Scalar:
        return "Scalar";
      case FieldInfo::TensorRank::Vector:
        return "Vector";
      default:
        return "Unknown TensorRank";
    }
}

// TODO: Fix style. This struct needn't be over 100 lines. Validation can occur in a
// separate function.
/**
 * @brief Structure to hold the variable attributes of a field. This includes things like
 * the name, equation type, whether it's nonlinear, and its dependence on other variables.
 */
struct VariableAttributes
{
  using EvalFlags = dealii::EvaluationFlags::EvaluationFlags;

  /**
   * @brief Allow VariableAttributeLoader to access private members.
   */
  friend class VariableAttributeLoader;

  FieldInfo field_info;

  /**
   * @brief Get the field name
   */
  [[nodiscard]] const std::string &
  get_name() const
  {
    return name;
  }

  /**
   * @brief Get the field index
   */
  [[nodiscard]] const Types::Index &
  get_field_index() const
  {
    return field_index;
  }

  /**
   * @brief Get the PDE type
   */
  [[nodiscard]] const PDEType &
  get_pde_type() const
  {
    return pde_type;
  }

  /**
   * @brief Whether the field is a postprocess variable.
   */
  [[nodiscard]] bool
  is_postprocess() const
  {
    return is_postprocessed_variable;
  }

  /**
   * @brief Whether the field is a nucleation rate variable.
   */
  [[nodiscard]] bool
  is_nucleation_rate() const
  {
    return is_nucleation_rate_variable;
  }

  /**
   * @brief Get the nucleating field indices.
   */
  [[nodiscard]] const std::vector<Types::Index> &
  get_nucleating_field_indices() const
  {
    return nucleation_indices;
  }

  /**
   * @brief Get the solve block.
   */
  [[nodiscard]] Types::Index
  get_solve_block() const
  {
    return solve_block;
  }

#ifdef ADDITIONAL_OPTIMIZATIONS
  /**
   * @brief Set the degenerate field index.
   */
  void
  set_degenerate_field_index(const Types::Index &_degenerate_field_index)
  {
    degenerate_field_index = _degenerate_field_index;
  }

  /**
   * @brief Get the degenerate field index.
   */
  [[nodiscard]] Types::Index
  get_degenerate_field_index() const
  {
    return degenerate_field_index;
  }
#endif

  /**
   * @brief Get the field solve type.
   */
  [[nodiscard]] FieldSolveType
  get_field_solve_type() const
  {
    return field_solve_type;
  }

  /**
   * @brief Get the RHS evaluation flags.
   */
  [[nodiscard]] const std::vector<std::vector<EvalFlags>> &
  get_eval_flag_set_rhs() const
  {
    return eval_flag_set_rhs;
  }

  /**
   * @brief Get the LHS evaluation flags.
   */
  [[nodiscard]] const std::vector<std::vector<EvalFlags>> &
  get_eval_flag_set_lhs() const
  {
    return eval_flag_set_lhs;
  }

  /**
   * @brief Get the RHS evaluation flags.
   */
  [[nodiscard]] const EvalFlags &
  get_eval_flags_residual_rhs() const
  {
    return eval_flags_residual_rhs;
  }

  /**
   * @brief Get the LHS evaluation flags.
   */
  [[nodiscard]] const EvalFlags &
  get_eval_flags_residual_lhs() const
  {
    return eval_flags_residual_lhs;
  }

  /**
   * @brief Get the dependency set for the RHS.
   */
  [[nodiscard]] const std::vector<std::vector<FieldInfo::TensorRank>> &
  get_dependency_set_rhs() const
  {
    return dependency_set_rhs;
  }

  /**
   * @brief Get the dependency set for the LHS.
   */
  [[nodiscard]] const std::vector<std::vector<FieldInfo::TensorRank>> &
  get_dependency_set_lhs() const
  {
    return dependency_set_lhs;
  }

  /**
   * @brief Get the max number of fields.
   */
  [[nodiscard]] Types::Index
  get_max_fields() const
  {
    Assert(max_fields != Numbers::invalid_index, dealii::ExcNotInitialized());
    return max_fields;
  }

  /**
   * @brief Get the max number of dependency types.
   */
  [[nodiscard]] Types::Index
  get_max_dependency_types() const
  {
    Assert(max_dependency_types != Numbers::invalid_index, dealii::ExcNotInitialized());
    return max_dependency_types;
  }

  /**
   * @brief Print variable attributes to summary.log
   */
  void
  print() const;

private:
  /**
   * @brief Combine 'value' and 'gradient' residual dependencies to one dependency set per
   * RHS and LHS. This will populate `dependencies_rhs` and `dependencies_lhs`.
   */
  void
  format_dependencies();

  /**
   * @brief Take user-defined dependency sets to set the residual flags for each
   * variable.
   */
  void
  parse_residual_dependencies();

  /**
   * @brief Take user-defined dependency sets to set the evaluation flags for each
   * variable.
   *
   * @param[in] other_var_attributes The other variable attributes.
   * @param[in] _max_fields The max number of fields that user has defined.
   * @param[in] _max_dependency_types The max number of dependency types.
   */
  void
  parse_dependencies(std::map<unsigned int, VariableAttributes> &other_var_attributes,
                     const Types::Index                         &_max_fields,
                     const Types::Index                         &_max_dependency_types);

  /**
   * @brief Using the assigned evaluation flags determine the solve type for this
   * equation.
   */
  void
  determine_field_solve_type(
    const std::map<unsigned int, VariableAttributes> &other_var_attributes);

  /**
   * @brief Validate a dependency.
   */
  void
  validate_dependency(const std::string           &variation,
                      DependencyType               dep_type,
                      const unsigned int          &other_index,
                      const FieldInfo::TensorRank &other_field_type,
                      const std::string           &context) const;

  /**
   * @brief Compute the dependency sets from eval_flag_set_rhs &
   * eval_flag_set_lhs
   */
  void
  compute_dependency_set(
    const std::map<unsigned int, VariableAttributes> &other_var_attributes);

  /**
   * @brief Compute the simplified dependency set from eval_flag_set_rhs &
   * eval_flag_set_lhs
   */
  void
  compute_simplified_dependency_set(
    const std::map<unsigned int, VariableAttributes> &other_var_attributes);

  /**
   * @brief Find the circular dependencies based on simple DFS algorithm.
   */
  void
  find_circular_dependencies(
    const std::map<unsigned int, VariableAttributes> &other_var_attributes);

  /**
   * @brief Recursive DFS
   */
  void
  recursive_depth_first_search(
    const std::map<unsigned int, VariableAttributes> &other_var_attributes,
    std::set<unsigned int>                           &visited,
    std::set<unsigned int>                           &current_stack,
    const unsigned int                               &vertex);

  /**
   * @brief Field name.
   * @remark User-set
   */
  std::string name;

  /**
   * @brief Field index.
   * @remark User-set
   */
  Types::Index field_index = Numbers::invalid_index;

  /**
   * @brief PDE type (Explicit/NONEXPLICIT).
   * @remark User-set
   */
  PDEType pde_type = Numbers::invalid_pde_type;

  /**
   * @brief Postprocess variable.
   * @remark User-set
   */
  bool is_postprocessed_variable = false;

  /**
   * @brief Is a nucleation rate
   * @remark User-set
   */
  bool is_nucleation_rate_variable = false;

  /**
   * @brief If this is a nucleation rate, the indices of the nucleating fields
   * @remark User-set
   */
  std::vector<Types::Index> nucleation_indices;

  /**
   * @brief Solve block
   * @remark User-set
   */
  Types::Index solve_block = 0;

#ifdef ADDITIONAL_OPTIMIZATIONS
  /**
   * @brief Degenerate field index.
   * @remark Internally determined
   */
  mutable Types::Index degenerate_field_index = Numbers::invalid_index;
#endif

  /**
   * @brief Internal classification for the field solve type.
   * @remark Internally determined
   */
  FieldSolveType field_solve_type = Numbers::invalid_field_solve_type;

  /**
   * @brief Max number of fields.
   */
  Types::Index max_fields = Numbers::invalid_index;

  /**
   * @brief Max number of dependency types.
   */
  Types::Index max_dependency_types = Numbers::invalid_index;

  /**
   * @brief A map of evaluation flags for the dependencies of the variable's RHS. This
   * will tell deal.II whether to evaluate the value, gradient, and/or hessian for the
   * specified field.
   *
   * The first dimension is the field index, the second dimension is the dependency type.
   * @remark Internally determined
   */
  std::vector<std::vector<EvalFlags>> eval_flag_set_rhs;

  /**
   * @brief A map of evaluation flags for the dependencies of the variable's LHS. This
   * will tell deal.II whether to evaluate the value, gradient, and/or hessian for the
   * specified field.
   *
   * The first dimension is the field index, the second dimension is the dependency type.
   * @remark Internally determined
   */
  std::vector<std::vector<EvalFlags>> eval_flag_set_lhs;

  /**
   * @brief Evaluation flags for the types of residual the user is expected to submit to
   * on the RHS.
   * @remark Internally determined
   */
  EvalFlags eval_flags_residual_rhs = dealii::EvaluationFlags::EvaluationFlags::nothing;

  /**
   * @brief Evaluation flags for the types of residual the user is expected to submit to
   * on the LHS. This is empty for Explicit fields.
   * @remark Internally determined
   */
  EvalFlags eval_flags_residual_lhs = dealii::EvaluationFlags::EvaluationFlags::nothing;

  /**
   * @brief TODO (Landinjm): Add brief description.
   *
   * To create the FEEvaluation object, we need to know two things: the index that
   * corresponds to the matrix-free object and the spacedim of the field (e.g., scalar or
   * vector). Notably, to avoid some burden of creating extra matrix-free objects for
   * multigrid, we occasionly use local indexing rather than a global index.
   *
   * The goal of this dependency set is to create a vector of FEEvaluation objects in
   * `VariableContainer`. For performance reasons, we represent this as a vector.
   *
   * TODO (Landinjm): Add some more details about the vector representation.
   *
   * @remark Internally determined
   */
  std::vector<std::vector<FieldInfo::TensorRank>> dependency_set_rhs;

  /**
   * @brief TODO (Landinjm): Add brief description.
   *
   * To create the FEEvaluation object, we need to know two things: the index that
   * corresponds to the matrix-free object and the spacedim of the field (e.g., scalar or
   * vector). Notably, to avoid some burden of creating extra matrix-free objects for
   * multigrid, we occasionly use local indexing rather than a global index.
   *
   * The goal of this dependency set is to create a vector of FEEvaluation objects in
   * `VariableContainer`. For performance reasons, we represent this as a vector.
   *
   * TODO (Landinjm): Add some more details about the vector representation.
   *
   * @remark Internally determined
   */
  std::vector<std::vector<FieldInfo::TensorRank>> dependency_set_lhs;

  /**
   * @brief Raw dependencies of the field.
   */
  RawDependencies raw_dependencies;

  /**
   * @brief A simplified set of evaluation flags for the dependencies of the variable's
   * LHS & RHS. This will help determine the FieldSolveType of the field. Fields indices
   * where the evaluation flags are not 0 (not nothing) are included. Additionally,
   * explicit fields are excluded to speed up the graph search.
   * @remark Internally determined
   */
  std::set<unsigned int> simplified_dependency_set;
};

PRISMS_PF_END_NAMESPACE
