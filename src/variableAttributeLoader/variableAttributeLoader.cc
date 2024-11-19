#include "../../include/variableAttributeLoader.h"

#include <deal.II/base/utilities.h>

#include <algorithm>

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