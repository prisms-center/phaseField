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
 * \brief Structure to hold the variable attributes of a field. This includes things like
 * the name, equation type, whether it's nonlinear, and its dependence on other variables.
 */
struct variableAttributes
{
  /**
   * \brief Field name. \remark User-set
   */
  std::string name;

  /**
   * \brief Field index. \remark User-set
   */
  types::index field_index = numbers::invalid_index;

  /**
   * \brief Field type (SCALAR/VECTOR). \remark User-set
   */
  fieldType field_type = fieldType::UNDEFINED_FIELD;

  /**
   * \brief PDE type (EXPLICIT/NONEXPLICIT). \remark User-set
   */
  PDEType pde_type = PDEType::UNDEFINED_PDE;

  /**
   * \brief Postprocess variable. \remark User-set
   */
  bool is_postprocess = false;

#ifdef ADDITIONAL_OPTIMIZATIONS
  /**
   * \brief Duplicate field index. \remark Internally determined
   */
  mutable types::index duplicate_field_index = numbers::invalid_index;
#endif

  /**
   * \brief Internal classification for the field solve type. \remark Internally
   * determined
   */
  fieldSolveType field_solve_type = fieldSolveType::UNDEFINED_SOLVE;

  /**
   * \brief The user-inputted dependencies for the RHS value term. \remark User-set
   */
  std::set<std::string> dependencies_value_RHS;

  /**
   * \brief The user-inputted dependencies for the RHS gradient term. \remark User-set
   */
  std::set<std::string> dependencies_gradient_RHS;

  /**
   * \brief The collection of value and gradient dependencies for the RHS. \remark
   * Internally determined
   */
  std::set<std::string> dependencies_RHS;

  /**
   * \brief The user-inputted dependencies for the LHS value term. \remark User-set
   */
  std::set<std::string> dependencies_value_LHS;

  /**
   * \brief The user-inputted dependencies for the LHS gradient term. \remark User-set
   */
  std::set<std::string> dependencies_gradient_LHS;

  /**
   * \brief The collection of value and gradient dependencies for the LHS. \remark
   * Internally determined
   */
  std::set<std::string> dependencies_LHS;

  /**
   * \brief A map of evaluation flags for the dependencies of the current variable's RHS.
   * This will tell deal.II whether to evaluate the value, gradient, and/or hessian for
   * the specified field. \remark Internally determined
   */
  std::map<std::pair<unsigned int, dependencyType>,
           dealii::EvaluationFlags::EvaluationFlags>
    eval_flag_set_RHS;

  /**
   * \brief A map of evaluation flags for the dependencies of the current variable's LHS.
   * This will tell deal.II whether to evaluate the value, gradient, and/or hessian for
   * the specified field. \remark Internally determined
   */
  std::map<std::pair<unsigned int, dependencyType>,
           dealii::EvaluationFlags::EvaluationFlags>
    eval_flag_set_LHS;

  /**
   * \brief Evaluation flags for the types of residual the user is expected to submit to
   * on the RHS. \remark Internally determined
   */
  dealii::EvaluationFlags::EvaluationFlags eval_flags_residual_RHS =
    dealii::EvaluationFlags::nothing;

  /**
   * \brief Evaluation flags for the types of residual the user is expected to submit to
   * on the LHS. This is empty for EXPLICIT fields. \remark Internally determined
   */
  dealii::EvaluationFlags::EvaluationFlags eval_flags_residual_LHS =
    dealii::EvaluationFlags::nothing;

  /**
   * \brief A dependency set where the RHS evaluation flags that are not 0 (not nothing)
   * are included. This is used to determine what FEEvaluatiob objects are necessary in
   * variable container. \remark Internally determined
   */
  std::map<unsigned int, std::map<dependencyType, fieldType>> dependency_set_RHS;

  /**
   * \brief A dependency set where the LHS evaluation flags that are not 0 (not nothing)
   * are included. This is used to determine what FEEvaluatiob objects are necessary in
   * variable container. \remark Internally determined
   */
  std::map<unsigned int, std::map<dependencyType, fieldType>> dependency_set_LHS;

  /**
   * \brief A simplified set of evaluation flags for the dependencies of the current
   * variable's LHS & RHS. This will help determine the fieldSolveType of the field.
   * Fields indices where the evaluation flags are not 0 (not nothing) are included.
   * Additionally, explicit fields are excluded to speed up the graph search. \remark
   * Internally determined
   */
  std::set<unsigned int> simplified_dependency_set;

  /**
   * \brief Combine 'value' and 'gradient' residual dependencies to one dependency set per
   * RHS and LHS. This will populate `dependencies_RHS` and `dependencies_LHS`.
   */
  void
  format_dependencies();

  /**
   * \brief Take user-defined dependency sets to set the residual flags for each
   * variable.
   */
  void
  parse_residual_dependencies();

  /**
   * \brief Take user-defined dependency sets to set the evaluation flags for each
   * variable.
   */
  void
  parse_dependencies(std::map<unsigned int, variableAttributes> &other_var_attributes);

  /**
   * \brief Using the assigned evaluation flags determine the solve type for this
   * equation.
   */
  void
  determine_field_solve_type(
    const std::map<unsigned int, variableAttributes> &other_var_attributes);

  /**
   * \brief Print variable attributes to summary.log
   */
  void
  print() const;

private:
  /**
   * \brief Validate a dependency.
   */
  void
  validate_dependency(const std::string  &variation,
                      dependencyType      dep_type,
                      const unsigned int &other_index,
                      const std::string  &context) const;

  /**
   * \brief Compute the dependency sets from eval_flag_set_RHS &
   * eval_flag_set_LHS
   */
  void
  compute_dependency_set(
    const std::map<unsigned int, variableAttributes> &other_var_attributes);

  /**
   * \brief Compute the simplified dependency set from eval_flag_set_RHS &
   * eval_flag_set_LHS
   */
  void
  compute_simplified_dependency_set(
    const std::map<unsigned int, variableAttributes> &other_var_attributes);

  /**
   * \brief Find the circular dependencies based on simple DFS algorithm.
   */
  void
  find_circular_dependencies(
    const std::map<unsigned int, variableAttributes> &other_var_attributes);

  /**
   * \brief Recursive DFS
   */
  void
  recursive_DFS(const std::map<unsigned int, variableAttributes> &other_var_attributes,
                std::set<unsigned int>                           &visited,
                std::set<unsigned int>                           &current_stack,
                const unsigned int                               &vertex);
};

PRISMS_PF_END_NAMESPACE
