#include "../../include/variableAttributeLoader.h"

#include <deal.II/base/mpi.h>

#include <algorithm>

void
variableAttributeLoader::format_dependencies()
{
  {
    for (auto &[index, variable] : attributes)
      {
        variable.dependencies_RHS.insert(variable.dependencies_value_RHS.begin(),
                                         variable.dependencies_value_RHS.end());
        variable.dependencies_RHS.insert(variable.dependencies_gradient_RHS.begin(),
                                         variable.dependencies_gradient_RHS.end());

        variable.dependencies_LHS.insert(variable.dependencies_value_LHS.begin(),
                                         variable.dependencies_value_LHS.end());
        variable.dependencies_LHS.insert(variable.dependencies_gradient_LHS.begin(),
                                         variable.dependencies_gradient_LHS.end());

        variable.dependencies_PP.insert(variable.dependencies_value_PP.begin(),
                                        variable.dependencies_value_PP.end());
        variable.dependencies_PP.insert(variable.dependencies_gradient_PP.begin(),
                                        variable.dependencies_gradient_PP.end());

        variable.dependency_set.insert(variable.dependencies_RHS.begin(),
                                       variable.dependencies_RHS.end());
        variable.dependency_set.insert(variable.dependencies_LHS.begin(),
                                       variable.dependencies_LHS.end());
        variable.dependency_set.insert(variable.dependencies_PP.begin(),
                                       variable.dependencies_PP.end());
      }
  }
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
                  for (auto &eval_flag : eval_flags_for_eq_type(other_variable.eq_type))
                    {
                      *eval_flag |= relevant_flag.at(variation);
                    }
                  other_variable.is_nonlinear |=
                    (eq_type != EXPLICIT_TIME_DEPENDENT) && (&other_variable == this) &&
                    (other_variable.eq_type != EXPLICIT_TIME_DEPENDENT);
                }
            }
        }
    }
}

// Constructor
variableAttributeLoader::variableAttributeLoader()
{
  setting_primary_field_attributes = true;
  loadVariableAttributes(); // This is the user-facing function
  setting_primary_field_attributes = false;
  loadPostProcessorVariableAttributes();

  format_dependencies();
  validate_attributes();
  for (auto &[index, variable] : attributes)
    {
      variable.parse_residual_dependencies();
      variable.parse_dependencies(attributes);
      variable.parse_dependencies(pp_attributes);
    }
  for (auto &[index, variable] : pp_attributes)
    {
      variable.parse_residual_dependencies();
    }
}

// Methods to set the various variable attributes
void
variableAttributeLoader::set_variable_name(const unsigned int &index,
                                           const std::string  &name)
{
  if (setting_primary_field_attributes)
    {
      attributes[index].name = name;
    }
  else
    {
      attributes[index].name_pp = name;
    }
}

void
variableAttributeLoader::set_variable_type(const unsigned int &index,
                                           const fieldType    &var_type)
{
  if (setting_primary_field_attributes)
    {
      attributes[index].var_type = var_type;
    }
  else
    {
      attributes[index].var_type_PP = var_type;
    }
}

void
variableAttributeLoader::set_variable_equation_type(const unsigned int &index,
                                                    const PDEType      &var_eq_type)
{
  attributes[index].eq_type = var_eq_type;
}

void
variableAttributeLoader::set_need_value_nucleation(const unsigned int &index,
                                                   const bool         &flag)
{
  attributes[index].need_value_nucleation = flag;
}

void
variableAttributeLoader::set_allowed_to_nucleate(const unsigned int &index,
                                                 const bool         &flag)
{
  attributes[index].nucleating_variable = flag;
}

void
variableAttributeLoader::set_output_integral(const unsigned int &index, const bool &flag)
{
  attributes[index].output_integral = flag;
}

void
variableAttributeLoader::set_dependencies_value_term_RHS(const unsigned int &index,
                                                         const std::string  &dependencies)
{
  std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  if (setting_primary_field_attributes)
    {
      attributes[index].dependencies_value_RHS =
        std::set<std::string>(dependencies_set.begin(), dependencies_set.end());
    }
  else
    {
      attributes[index].dependencies_value_PP =
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
  if (setting_primary_field_attributes)
    {
      attributes[index].dependencies_gradient_RHS =
        std::set<std::string>(dependencies_set.begin(), dependencies_set.end());
    }
  else
    {
      attributes[index].dependencies_gradient_PP =
        std::set<std::string>(dependencies_set.begin(), dependencies_set.end());
    }
}

void
variableAttributeLoader::set_dependencies_value_term_LHS(const unsigned int &index,
                                                         const std::string  &dependencies)
{
  std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  attributes[index].dependencies_value_LHS =
    std::set<std::string>(dependencies_set.begin(), dependencies_set.end());
}

void
variableAttributeLoader::set_dependencies_gradient_term_LHS(
  const unsigned int &index,
  const std::string  &dependencies)
{
  std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  attributes[index].dependencies_gradient_LHS =
    std::set<std::string>(dependencies_set.begin(), dependencies_set.end());
}

void
variableAttributeLoader::validate_attributes()
{
  for (const auto &[index, attribute_set] : attributes)
    {
    }
}

std::string
variableAttributeLoader::strip_whitespace(const std::string &_text)
{
  std::string text = _text;
  text.erase(std::remove(text.begin(), text.end(), ' '), text.end());
  return text;
}