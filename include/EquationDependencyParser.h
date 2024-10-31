#ifndef INCLUDE_EQUATIONDEPENDECYPARSER_H_
#define INCLUDE_EQUATIONDEPENDECYPARSER_H_

#include <deal.II/matrix_free/evaluation_flags.h>

#include "varTypeEnums.h"

#include <algorithm>
#include <string>
#include <vector>

/**
 * This is a class to parse the strings the user puts into the
 * variableAttributeLoader to specify which variable values, gradients,
 * hessians, etc are needed for each governing equation.
 */
class EquationDependencyParser
{
public:
  void
  parse(std::vector<std::string> &var_name,
        std::vector<PDEType>      var_eq_type,
        std::vector<std::string>  sorted_dependencies_value_RHS,
        std::vector<std::string>  sorted_dependencies_gradient_RHS,
        std::vector<std::string>  sorted_dependencies_value_LHS,
        std::vector<std::string>  sorted_dependencies_gradient_LHS,
        std::vector<bool>        &var_nonlinear);

  void
  pp_parse(std::vector<std::string> &var_name,
           std::vector<std::string> &pp_var_name,
           std::vector<std::string>  sorted_dependencies_value,
           std::vector<std::string>  sorted_dependencies_gradient);

  // Evaluation flags for each type of solution variable (e.g., explicit, nonexplicit,
  // nonexplicit_change, etc.)
  std::vector<dealii::EvaluationFlags::EvaluationFlags> eval_flags_explicit_RHS;
  std::vector<dealii::EvaluationFlags::EvaluationFlags> eval_flags_nonexplicit_RHS;
  std::vector<dealii::EvaluationFlags::EvaluationFlags> eval_flags_nonexplicit_LHS;
  std::vector<dealii::EvaluationFlags::EvaluationFlags> eval_flags_change_nonexplicit_LHS;

  std::vector<dealii::EvaluationFlags::EvaluationFlags> eval_flags_residual_explicit_RHS;
  std::vector<dealii::EvaluationFlags::EvaluationFlags>
    eval_flags_residual_nonexplicit_RHS;
  std::vector<dealii::EvaluationFlags::EvaluationFlags>
    eval_flags_residual_nonexplicit_LHS;

  // All of the vectors of flags for what is needed for the postprocessing
  // variables
  std::vector<dealii::EvaluationFlags::EvaluationFlags> eval_flags_postprocess;
  std::vector<dealii::EvaluationFlags::EvaluationFlags> eval_flags_residual_postprocess;

protected:
  /*
   * Method to strip excess whitespace for the dependency lists
   */
  void
  strip_dependency_whitespace(std::string &dependency_list);

  /**
   * Method to parse the RHS dependency strings and populate the vectors for
   * whether values, gradients, or hessians are needed.
   */
  void
  parseDependencyListRHS(
    std::vector<std::string>                              &variable_name_list,
    std::vector<PDEType>                                   variable_eq_type,
    unsigned int                                           variable_index,
    std::string                                           &value_dependencies,
    std::string                                           &gradient_dependencies,
    std::vector<dealii::EvaluationFlags::EvaluationFlags> &evaluation_flags,
    std::vector<dealii::EvaluationFlags::EvaluationFlags> &residual_flags,
    bool                                                  &is_nonlinear);

  /**
   * Method to parse the LHS dependency strings and populate the vectors for
   * whether values, gradients, or hessians are needed.
   */
  void
  parseDependencyListLHS(
    std::vector<std::string>                              &variable_name_list,
    std::vector<PDEType>                                   variable_eq_type,
    unsigned int                                           variable_index,
    std::string                                           &value_dependencies,
    std::string                                           &gradient_dependencies,
    std::vector<dealii::EvaluationFlags::EvaluationFlags> &evaluation_flags,
    std::vector<dealii::EvaluationFlags::EvaluationFlags> &change_flags,
    std::vector<dealii::EvaluationFlags::EvaluationFlags> &residual_flags,
    bool                                                  &is_nonlinear);

  /**
   * Method to parse the postprocessing dependency strings and populate the
   * vectors for whether values, gradients, or hessians are needed.
   */
  void
  parseDependencyListPP(
    std::vector<std::string>                              &variable_name_list,
    unsigned int                                           variable_index,
    std::string                                           &value_dependencies,
    std::string                                           &gradient_dependencies,
    std::vector<dealii::EvaluationFlags::EvaluationFlags> &evaluation_flags,
    std::vector<dealii::EvaluationFlags::EvaluationFlags> &residual_flags);
};

#endif
