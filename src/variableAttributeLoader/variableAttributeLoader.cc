#include "../../include/variableAttributeLoader.h"

#include <deal.II/base/mpi.h>

#include "varTypeEnums.h"

#include <algorithm>

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

// Constructor
variableAttributeLoader::variableAttributeLoader()
{
  relevant_attributes = &attributes;
  loadVariableAttributes(); // This is a user-facing function
  relevant_attributes = &pp_attributes;
  loadPostProcessorVariableAttributes(); // This is a user-facing function
  relevant_attributes = nullptr;

  for (auto &[index, pp_variable] : pp_attributes)
    {
      pp_variable.is_pp   = true;
      pp_variable.eq_type = EXPLICIT_TIME_DEPENDENT;
    }

  for (auto &[index, variable] : attributes)
    {
      variable.format_dependencies();
    }
  for (auto &[pp_index, pp_variable] : pp_attributes)
    {
      pp_variable.format_dependencies();
    }
  validate_attributes();
  for (auto &[index, variable] : attributes)
    {
      variable.parse_residual_dependencies();
      variable.parse_dependencies(attributes);
      variable.parse_dependencies(pp_attributes);
    }
  for (auto &[index, pp_variable] : pp_attributes)
    {
      pp_variable.parse_residual_dependencies();
    }
}

// Methods to set the various variable attributes
void
variableAttributeLoader::set_variable_name(const unsigned int &index,
                                           const std::string  &name)
{
  (*relevant_attributes)[index].name = name;
}

void
variableAttributeLoader::set_variable_type(const unsigned int &index,
                                           const fieldType    &var_type)
{
  (*relevant_attributes)[index].var_type = var_type;
}

void
variableAttributeLoader::set_variable_equation_type(const unsigned int &index,
                                                    const PDEType      &var_eq_type)
{
  (*relevant_attributes)[index].eq_type = var_eq_type;
}

void
variableAttributeLoader::set_need_value_nucleation(const unsigned int &index,
                                                   const bool         &flag)
{
  (*relevant_attributes)[index].need_value_nucleation = flag;
}

void
variableAttributeLoader::set_allowed_to_nucleate(const unsigned int &index,
                                                 const bool         &flag)
{
  (*relevant_attributes)[index].nucleating_variable = flag;
}

void
variableAttributeLoader::set_output_integral(const unsigned int &index, const bool &flag)
{
  (*relevant_attributes)[index].output_integral = flag;
  (*relevant_attributes)[index].calc_integral   = flag;
}

void
variableAttributeLoader::set_dependencies_value_term_RHS(const unsigned int &index,
                                                         const std::string  &dependencies)
{
  std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  /* (*relevant_attributes)[index].dependencies_value_RHS =
      std::set<std::string>(dependencies_set.begin(), dependencies_set.end()); */
  if (relevant_attributes != &pp_attributes)
    {
      attributes[index].dependencies_value_RHS =
        std::set<std::string>(dependencies_set.begin(), dependencies_set.end());
    }
  else
    {
      pp_attributes[index].dependencies_value_PP =
        std::set<std::string>(dependencies_set.begin(), dependencies_set.end());
    }
}

void
variableAttributeLoader::set_dependencies_gradient_term_RHS(
  const unsigned int &index,
  const std::string  &dependencies)
{
  std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  /* (*relevant_attributes)[index].dependencies_gradient_RHS =
    std::set<std::string>(dependencies_set.begin(), dependencies_set.end()); */
  if (relevant_attributes != &pp_attributes)
    {
      attributes[index].dependencies_gradient_RHS =
        std::set<std::string>(dependencies_set.begin(), dependencies_set.end());
    }
  else
    {
      pp_attributes[index].dependencies_gradient_PP =
        std::set<std::string>(dependencies_set.begin(), dependencies_set.end());
    }
}

void
variableAttributeLoader::set_dependencies_value_term_LHS(const unsigned int &index,
                                                         const std::string  &dependencies)
{
  std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  (*relevant_attributes)[index].dependencies_value_LHS =
    std::set<std::string>(dependencies_set.begin(), dependencies_set.end());
}

void
variableAttributeLoader::set_dependencies_gradient_term_LHS(
  const unsigned int &index,
  const std::string  &dependencies)
{
  std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  (*relevant_attributes)[index].dependencies_gradient_LHS =
    std::set<std::string>(dependencies_set.begin(), dependencies_set.end());
}

void
variableAttributeLoader::validate_attributes()
{
  // Make sure main attributes arent set in pp_attributes
  // Make sure dependencies are all variable names and if there are "change()" deps, they
  // are correctly placed
  /* for (const auto &[index, attribute_set] : attributes)
    {
    }
  for (const auto &[pp_index, pp_attribute_set] : pp_attributes)
    {
    } */
}

std::string
variableAttributeLoader::strip_whitespace(const std::string &_text)
{
  std::string text = _text;
  text.erase(std::remove(text.begin(), text.end(), ' '), text.end());
  return text;
}