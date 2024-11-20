// Class to hold the variable attributes that will be passed to a
// userInputParameters object
#ifndef VARIABLEATTRIBUTES_H
#define VARIABLEATTRIBUTES_H

#include <deal.II/matrix_free/evaluation_flags.h>

#include "varTypeEnums.h"

#include <map>
#include <set>
#include <string>

using EvalFlags = dealii::EvaluationFlags::EvaluationFlags;

struct variableAttributes
{
  // Variable attributes
  std::string name                  = "";
  fieldType   var_type              = UNDEFINED_FIELD;
  PDEType     eq_type               = UNDEFINED_PDE;
  bool        need_value_nucleation = false;
  bool        nucleating_variable   = false;
  bool        is_pp                 = false;
  bool        is_nonlinear          = false;
  bool        calc_integral         = false;
  bool        output_integral       = false;

  // This variable's dependencies
  std::set<std::string> dependencies_value_RHS;
  std::set<std::string> dependencies_gradient_RHS;
  std::set<std::string> dependencies_RHS;
  std::set<std::string> dependencies_value_LHS;
  std::set<std::string> dependencies_gradient_LHS;
  std::set<std::string> dependencies_LHS;
  std::set<std::string> dependencies_value_PP;
  std::set<std::string> dependencies_gradient_PP;
  std::set<std::string> dependencies_PP;

  std::set<std::string> dependency_set;

  // Evaluation Flags. Tells deal.ii whether or not to retrieve the value, grad, hess,
  // etc. for equations
  EvalFlags eval_flags_explicit_RHS    = dealii::EvaluationFlags::nothing;
  EvalFlags eval_flags_nonexplicit_RHS = dealii::EvaluationFlags::nothing;
  EvalFlags eval_flags_nonexplicit_LHS = dealii::EvaluationFlags::nothing;

  EvalFlags eval_flags_change_nonexplicit_LHS = dealii::EvaluationFlags::nothing;

  EvalFlags eval_flags_residual_explicit_RHS    = dealii::EvaluationFlags::nothing;
  EvalFlags eval_flags_residual_nonexplicit_RHS = dealii::EvaluationFlags::nothing;
  EvalFlags eval_flags_residual_nonexplicit_LHS = dealii::EvaluationFlags::nothing;

  EvalFlags eval_flags_postprocess          = dealii::EvaluationFlags::nothing;
  EvalFlags eval_flags_residual_postprocess = dealii::EvaluationFlags::nothing;

  /**
   * \brief Combine 'value' and 'gradient' dependencies to one dependency set per RHS,
   * LHS, PP.
   */
  void
  format_dependencies();

  /**
   * \brief Take user-defined dependency sets to set the evaluation flags for each
   * variable.
   */
  void
  parse_dependencies(std::map<uint, variableAttributes> &other_var_attributes);

  /**
   * \brief Take user-defined dependency sets to set the residual flags for each
   * variable.
   */
  void
  parse_residual_dependencies();

  /**
   * \brief Helper function that returns a set of pointers to the flags that need to be
   * set when `other_variable` is dependent on this variable.
   * \param other_variable Variable that is dependent on this variable.
   */
  std::set<EvalFlags *>
  eval_flags_for_eq_type(const variableAttributes &other_variable);

  /*   std::set<EvalFlags *>
    residual_flags_for_eq_type(PDEType other_eq_type)
    {
      if (other_eq_type == EXPLICIT_TIME_DEPENDENT)
        {
          return {&eval_flags_residual_explicit_RHS};
        }
      if (other_eq_type == AUXILIARY)
        {
          return {&eval_flags_residual_nonexplicit_RHS};
        }
      if (other_eq_type == IMPLICIT_TIME_DEPENDENT || other_eq_type == TIME_INDEPENDENT)
        {
          return {&eval_flags_residual_nonexplicit_RHS,
                  &eval_flags_residual_nonexplicit_LHS};
        }
      return {};
    } */

  /* std::set<std::pair<std::set<std::string> *, std::set<std::string> *>>
    dep_set_for_eq_type(PDEType other_eq_type)
    {
      if (other_eq_type == EXPLICIT_TIME_DEPENDENT)
        {
          return dependencies_RHS;
          // return {
          //   {&dependencies_value_RHS, &dependencies_gradient_RHS}
          // };
        }
      if (other_eq_type == AUXILIARY)
        {
          return dependencies_RHS;
          // return {
          //   {&dependencies_value_RHS, &dependencies_gradient_RHS}
          // };
        }
      if (other_eq_type == IMPLICIT_TIME_DEPENDENT || other_eq_type == TIME_INDEPENDENT)
        {
          return dependencies_RHS;
          // return {
          //   {&dependencies_value_RHS, &dependencies_gradient_RHS}
          // };
        }
      return {};
    } */
};

#endif
