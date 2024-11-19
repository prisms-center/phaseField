#include "../../include/variableAttributes.h"

void
variableAttributes::format_dependencies()
{
  dependencies_RHS.insert(dependencies_value_RHS.begin(), dependencies_value_RHS.end());
  dependencies_RHS.insert(dependencies_gradient_RHS.begin(),
                          dependencies_gradient_RHS.end());

  dependencies_LHS.insert(dependencies_value_LHS.begin(), dependencies_value_LHS.end());
  dependencies_LHS.insert(dependencies_gradient_LHS.begin(),
                          dependencies_gradient_LHS.end());

  dependencies_PP.insert(dependencies_value_PP.begin(), dependencies_value_PP.end());
  dependencies_PP.insert(dependencies_gradient_PP.begin(),
                         dependencies_gradient_PP.end());

  dependency_set.insert(dependencies_RHS.begin(), dependencies_RHS.end());
  dependency_set.insert(dependencies_LHS.begin(), dependencies_LHS.end());
  dependency_set.insert(dependencies_PP.begin(), dependencies_PP.end());
}

void
variableAttributes::parse_residual_dependencies()
{ // Check if either is empty and set value and gradient flags for the
  // residual appropriately
  const std::map<
    PDEType,
    std::set<std::pair<EvalFlags *,
                       std::pair<std::set<std::string> *, std::set<std::string> *>>>>
    dependencies_for = {
      {EXPLICIT_TIME_DEPENDENT,
       {{&eval_flags_residual_explicit_RHS,
         {&dependencies_value_RHS, &dependencies_gradient_RHS}}}},
      {AUXILIARY,
       {{&eval_flags_residual_nonexplicit_RHS,
         {&dependencies_value_RHS, &dependencies_gradient_RHS}}}},
      {IMPLICIT_TIME_DEPENDENT,
       {{&eval_flags_residual_nonexplicit_RHS,
         {&dependencies_value_RHS, &dependencies_gradient_RHS}},
        {&eval_flags_residual_nonexplicit_LHS,
         {&dependencies_value_LHS, &dependencies_gradient_LHS}}}},
      {TIME_INDEPENDENT,
       {{&eval_flags_residual_nonexplicit_RHS,
         {&dependencies_value_RHS, &dependencies_gradient_RHS}},
        {&eval_flags_residual_nonexplicit_LHS,
         {&dependencies_value_LHS, &dependencies_gradient_LHS}}}}
  };

  if (is_pp)
    {
      EvalFlags             &residual_flags          = eval_flags_residual_postprocess;
      std::set<std::string> &value_dependency_set    = dependencies_value_PP;
      std::set<std::string> &gradient_dependency_set = dependencies_gradient_PP;

      if (!value_dependency_set.empty())
        {
          residual_flags |= dealii::EvaluationFlags::values;
        }
      if (!gradient_dependency_set.empty())
        {
          residual_flags |= dealii::EvaluationFlags::gradients;
        }
      return;
    }

  for (const std::pair<EvalFlags *,
                       std::pair<std::set<std::string> *, std::set<std::string> *>>
         &eval_dep_pair : dependencies_for.at(eq_type))
    {
      EvalFlags             &residual_flags          = *(eval_dep_pair.first);
      auto                  &dep_pair                = eval_dep_pair.second;
      std::set<std::string> &value_dependency_set    = *(dep_pair.first);
      std::set<std::string> &gradient_dependency_set = *(dep_pair.second);

      if (!value_dependency_set.empty())
        {
          residual_flags |= dealii::EvaluationFlags::values;
        }
      if (!gradient_dependency_set.empty())
        {
          residual_flags |= dealii::EvaluationFlags::gradients;
        }
    }
}

void
variableAttributes::parse_dependencies(
  std::map<uint, variableAttributes> &other_var_attributes)
{
  const std::map<std::string, EvalFlags> relevant_flag {
    {"val",         dealii::EvaluationFlags::values   },
    {"grad",        dealii::EvaluationFlags::gradients},
    {"hess",        dealii::EvaluationFlags::hessians },
    {"val_change",  dealii::EvaluationFlags::values   },
    {"grad_change", dealii::EvaluationFlags::gradients},
    {"hess_change", dealii::EvaluationFlags::hessians }
  };
  const std::map<std::string, std::pair<std::string, std::string>> delimiters = {
    {"val",         {"", ""}              },
    {"grad",        {"grad(", ")"}        },
    {"hess",        {"hess(", ")"}        },
    {"val_change",  {"change(", ")"}      },
    {"grad_change", {"grad(change(", "))"}},
    {"hess_change", {"hess(change(", "))"}}
  };

  for (const auto &[variation, delimiter] : delimiters)
    {
      bool is_change = variation == "val_change" || variation == "grad_change" ||
                       variation == "hess_change";
      std::string possible_dependency = delimiter.first + name + delimiter.second;
      if (is_change && dependency_set.find(possible_dependency) != dependency_set.end())
        {
          eval_flags_change_nonexplicit_LHS |= relevant_flag.at(variation);
        }
      else
        {
          for (auto &[other_index, other_variable] : other_var_attributes)
            {
              if (other_variable.dependency_set.find(possible_dependency) !=
                  other_variable.dependency_set.end())
                {
                  for (auto &eval_flag : eval_flags_for_eq_type(other_variable))
                    {
                      *eval_flag |= relevant_flag.at(variation);
                    }
                  other_variable.is_nonlinear |=
                    (eq_type != EXPLICIT_TIME_DEPENDENT) &&
                    (other_variable.eq_type != EXPLICIT_TIME_DEPENDENT) &&
                    (&other_variable != this);
                }
            }
        }
    }
}

std::set<EvalFlags *>
variableAttributes::eval_flags_for_eq_type(const variableAttributes &other_variable)
{
  PDEType other_eq_type = other_variable.eq_type;
  if (other_variable.is_pp)
    {
      return {&eval_flags_postprocess};
    }
  if (other_eq_type == EXPLICIT_TIME_DEPENDENT)
    {
      return {&eval_flags_explicit_RHS};
    }
  if (other_eq_type == AUXILIARY)
    {
      return {&eval_flags_nonexplicit_RHS};
    }
  if (other_eq_type == IMPLICIT_TIME_DEPENDENT || other_eq_type == TIME_INDEPENDENT)
    {
      return {&eval_flags_nonexplicit_RHS, &eval_flags_nonexplicit_LHS};
    }
  return {};
}
