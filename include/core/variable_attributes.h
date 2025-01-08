#ifndef variable_attributes_h
#define variable_attributes_h

#include <deal.II/matrix_free/evaluation_flags.h>

#include <core/type_enums.h>
#include <map>
#include <set>
#include <string>
#include <unordered_map>

using EvalFlags = dealii::EvaluationFlags::EvaluationFlags;
struct variableAttributes;
using AttributesList = std::map<uint, variableAttributes>;

/**
 * \brief Simple hash function for pairs.
 */
struct pairHash
{
public:
  template <typename T, typename U>
  std::size_t
  operator()(const std::pair<T, U> &x) const
  {
    return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
  }
};

/**
 * \brief Structure to hold the variable attributes of a field. This includes things like
 * the name, equation type, whether it's nonlinear, and its dependence on other variables.
 */
struct variableAttributes
{
  /**
   * \brief Field name.
   */
  std::string name;

  /**
   * \brief Field index.
   */
  uint field_index;

  /**
   * \brief Field type (SCALAR/VECTOR).
   */
  fieldType field_type = fieldType::UNDEFINED_FIELD;

  /**
   * \brief PDE type (EXPLICIT/NONEXPLICIT).
   */
  PDEType pde_type = PDEType::UNDEFINED_PDE;

  /**
   * \brief Internal classification for the field solve type.
   */
  fieldSolveType field_solve_type = fieldSolveType::UNDEFINED_SOLVE;

  /**
   * \brief The user-inputted dependencies for the RHS value term.
   */
  std::set<std::string> dependencies_value_RHS;

  /**
   * \brief The user-inputted dependencies for the RHS gradient term.
   */
  std::set<std::string> dependencies_gradient_RHS;

  /**
   * \brief The collection of value and gradient dependencies for the RHS.
   */
  std::set<std::string> dependencies_RHS;

  /**
   * \brief The user-inputted dependencies for the LHS value term.
   */
  std::set<std::string> dependencies_value_LHS;

  /**
   * \brief The user-inputted dependencies for the LHS gradient term.
   */
  std::set<std::string> dependencies_gradient_LHS;

  /**
   * \brief The collection of value and gradient dependencies for the LHS.
   */
  std::set<std::string> dependencies_LHS;

  /**
   * \brief A map of evaluation flags for the dependencies of the current variable's RHS.
   * This will tell deal.II whether to evaluate the value, gradient, and/or hessian for
   * the specified field.
   */
  std::unordered_map<std::pair<uint, dependencyType>, EvalFlags, pairHash>
    eval_flag_set_RHS;

  /**
   * \brief A map of evaluation flags for the dependencies of the current variable's LHS.
   * This will tell deal.II whether to evaluate the value, gradient, and/or hessian for
   * the specified field.
   */
  std::unordered_map<std::pair<uint, dependencyType>, EvalFlags, pairHash>
    eval_flag_set_LHS;

  /**
   * \brief Evaluation flags for the types of residual the user is expected to submit to
   * on the RHS.
   */
  EvalFlags eval_flags_residual_RHS = dealii::EvaluationFlags::nothing;

  /**
   * \brief Evaluation flags for the types of residual the user is expected to submit to
   * on the LHS. This is empty for EXPLICIT fields.
   */
  EvalFlags eval_flags_residual_LHS = dealii::EvaluationFlags::nothing;

  /**
   * \brief A simplified set of evaluation flags for the dependencies of the current
   * variable's LHS & RHS. This will help determine the fieldSolveType of the field.
   * Fields indices where the evaluation flags are not 0 (not nothing) are included.
   * Additionally, explicit fields are excluded to speed up the graph search.
   */
  std::set<uint> dependency_set;

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
  parse_dependencies(AttributesList &other_var_attributes);

  /**
   * \brief Using the assigned evaluation flags determine the solve type for this
   * equation.
   */
  void
  determine_field_solve_type(AttributesList &other_var_attributes,
                             const bool     &is_postprocessed = false);

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
  validate_dependency(const std::string &variation,
                      dependencyType     dep_type,
                      const uint        &other_index,
                      const std::string &context) const;

  /**
   * \brief Compute the dependency set from eval_flag_set_RHS & eval_flag_set_LHS
   */
  void
  compute_dependency_set(const AttributesList &other_var_attributes);

  /**
   * \brief Find the circular dependencies based on simple DFS algorithm.
   */
  void
  find_circular_dependencies(const AttributesList &other_var_attributes);

  /**
   * \brief Recursive DFS
   */
  void
  recursive_DFS(const AttributesList &other_var_attributes,
                std::set<uint>       &visited,
                std::set<uint>       &current_stack,
                const uint           &vertex);
};

#endif
