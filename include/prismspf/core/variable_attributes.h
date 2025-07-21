// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include <map>
#include <set>
#include <string>
#include <utility>

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
};

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
   * @brief Get the field type
   */
  [[nodiscard]] const FieldType &
  get_field_type() const
  {
    return field_type;
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

#ifdef ADDITIONAL_OPTIMIZATIONS
  /**
   * @brief Set the duplicate field index.
   */
  void
  set_duplicate_field_index(const Types::Index &_duplicate_field_index)
  {
    duplicate_field_index = _duplicate_field_index;
  }

  /**
   * @brief Get the duplicate field index.
   */
  [[nodiscard]] Types::Index
  get_duplicate_field_index() const
  {
    return duplicate_field_index;
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
   * @brief Get the RHS evaluation flags.
   */
  [[nodiscard]] std::vector<std::vector<EvalFlags>> &
  get_eval_flag_set_rhs()
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
  [[nodiscard]] const std::vector<std::vector<FieldType>> &
  get_dependency_set_rhs() const
  {
    return dependency_set_rhs;
  }

  /**
   * @brief Get the dependency set for the RHS.
   */
  [[nodiscard]] std::vector<std::vector<FieldType>> &
  get_dependency_set_rhs()
  {
    return dependency_set_rhs;
  }

  /**
   * @brief Set the dependency set for the RHS.
   */
  void
  set_dependency_set_rhs(const std::vector<std::vector<FieldType>> &_dependency_set_rhs)
  {
    dependency_set_rhs = _dependency_set_rhs;
  }

  /**
   * @brief Get the dependency set for the LHS.
   */
  [[nodiscard]] const std::vector<std::vector<FieldType>> &
  get_dependency_set_lhs() const
  {
    return dependency_set_lhs;
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
   * @param[in] max_fields The max number of fields that user has defined.
   * @param[in] max_dependency_types The max number of dependency types.
   */
  void
  parse_dependencies(std::map<unsigned int, VariableAttributes> &other_var_attributes,
                     const Types::Index                         &max_fields,
                     const Types::Index                         &max_dependency_types);

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
  validate_dependency(const std::string  &variation,
                      DependencyType      dep_type,
                      const unsigned int &other_index,
                      const std::string  &context) const;

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
   * @brief Field type (Scalar/Vector).
   * @remark User-set
   */
  FieldType field_type = Numbers::invalid_field_type;

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

#ifdef ADDITIONAL_OPTIMIZATIONS
  /**
   * @brief Duplicate field index.
   * @remark Internally determined
   *
   * TODO (landinjm): Rename
   */
  mutable Types::Index duplicate_field_index = Numbers::invalid_index;
#endif

  /**
   * @brief Internal classification for the field solve type.
   * @remark Internally determined
   */
  FieldSolveType field_solve_type = Numbers::invalid_field_solve_type;

  /**
   * @brief A map of evaluation flags for the dependencies of the current variable's RHS.
   * This will tell deal.II whether to evaluate the value, gradient, and/or hessian for
   * the specified field.
   *
   * The first dimension is the field index, the second dimension is the dependency type.
   * @remark Internally determined
   */
  std::vector<std::vector<EvalFlags>> eval_flag_set_rhs;

  /**
   * @brief A map of evaluation flags for the dependencies of the current variable's LHS.
   * This will tell deal.II whether to evaluate the value, gradient, and/or hessian for
   * the specified field.
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
  std::vector<std::vector<FieldType>> dependency_set_rhs;

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
  std::vector<std::vector<FieldType>> dependency_set_lhs;

  /**
   * @brief Raw dependencies of the field.
   */
  RawDependencies raw_dependencies;

  /**
   * @brief A simplified set of evaluation flags for the dependencies of the current
   * variable's LHS & RHS. This will help determine the FieldSolveType of the field.
   * Fields indices where the evaluation flags are not 0 (not nothing) are included.
   * Additionally, explicit fields are excluded to speed up the graph search.
   * @remark Internally determined
   */
  std::set<unsigned int> simplified_dependency_set;
};

PRISMS_PF_END_NAMESPACE
