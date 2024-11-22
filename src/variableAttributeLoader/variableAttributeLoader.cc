#include "../../include/variableAttributeLoader.h"

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>

// NOLINTBEGIN(cppcoreguidelines-prefer-member-initializer)
variableAttributeLoader::variableAttributeLoader()
{
  relevant_attributes = &attributes;
  loadVariableAttributes();
  relevant_attributes = &pp_attributes;
  loadPostProcessorVariableAttributes();
  relevant_attributes = nullptr;

  for (auto &[pp_index, pp_variable] : pp_attributes)
    {
      pp_variable.is_pp = true;
      Assert(pp_variable.eq_type == EXPLICIT_TIME_DEPENDENT ||
               pp_variable.eq_type == UNDEFINED_PDE,
             dealii::ExcMessage("PRISMS-PF Error: Warning: Postprocess variables can "
                                "only be explicit.\nProblem postprocessing index: " +
                                std::to_string(pp_index) +
                                ", Variable name: " + pp_variable.name + "\n"));
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

// NOLINTEND(cppcoreguidelines-prefer-member-initializer)

// Methods to set the various variable attributes
void
variableAttributeLoader::set_variable_name(const unsigned int &index,
                                           const std::string  &name) const
{
  (*relevant_attributes)[index].name = name;
}

void
variableAttributeLoader::set_variable_type(const unsigned int &index,
                                           const fieldType    &var_type) const
{
  (*relevant_attributes)[index].var_type = var_type;
}

void
variableAttributeLoader::set_variable_equation_type(const unsigned int &index,
                                                    const PDEType      &var_eq_type) const
{
  (*relevant_attributes)[index].eq_type = var_eq_type;
}

void
variableAttributeLoader::set_need_value_nucleation(const unsigned int &index,
                                                   const bool         &flag) const
{
  (*relevant_attributes)[index].need_value_nucleation = flag;
}

void
variableAttributeLoader::set_allowed_to_nucleate(const unsigned int &index,
                                                 const bool         &flag) const
{
  (*relevant_attributes)[index].nucleating_variable = flag;
}

void
variableAttributeLoader::set_output_integral(const unsigned int &index,
                                             const bool         &flag) const
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
variableAttributeLoader::set_dependencies_value_term_LHS(
  const unsigned int &index,
  const std::string  &dependencies) const
{
  std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  (*relevant_attributes)[index].dependencies_value_LHS =
    std::set<std::string>(dependencies_set.begin(), dependencies_set.end());
}

void
variableAttributeLoader::set_dependencies_gradient_term_LHS(
  const unsigned int &index,
  const std::string  &dependencies) const
{
  std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  (*relevant_attributes)[index].dependencies_gradient_LHS =
    std::set<std::string>(dependencies_set.begin(), dependencies_set.end());
}

void
variableAttributeLoader::validate_variable_name(
  const std::string           &name,
  const std::set<std::string> &forbidden_names,
  const std::string           &context,
  unsigned int                 index)
{
  Assert(!name.empty(),
         dealii::ExcMessage("PRISMS-PF Error: " + context +
                            " Variable names must not be empty.\nProblem index: " +
                            std::to_string(index)));

  for (const std::string &forbidden_name : forbidden_names)
    {
      std::string error_message = "PRISMS-PF Error: " + context +
                                  " Variable names must not contain \"grad()\", "
                                  "\"hess()\", \"change()\".\nProblem index: ";
      error_message.append(std::to_string(index).append(", Variable name: " + name));

      Assert(name.find(forbidden_name) == std::string::npos,
             dealii::ExcMessage(error_message));
    }
}

void
variableAttributeLoader::populate_dependencies(
  const std::set<std::pair<std::string, std::string>> &reg_delimiters,
  const std::string                                   &variable_name,
  unsigned int                                         index,
  std::set<std::string>                               &reg_possible_deps,
  std::map<uint, std::set<std::string>>               &change_possible_deps)
{
  for (const auto &delims : reg_delimiters)
    {
      reg_possible_deps.insert(delims.first + variable_name + delims.second);
      change_possible_deps[index].insert(delims.first + "change(" + variable_name + ")" +
                                         delims.second);
    }
}

void
variableAttributeLoader::validate_dependencies(
  const std::set<std::string>                 &dependencies,
  const std::string                           &context,
  unsigned int                                 index,
  const std::string                           &variable_name,
  const std::set<std::string>                 &reg_possible_deps,
  const std::map<uint, std::set<std::string>> &change_possible_deps)
{
  for (const std::string &dependency : dependencies)
    {
      std::string error_message =
        "PRISMS-PF Error: Invalid " + context + " dependency.\nProblem index: ";
      error_message.append(
        std::to_string(index).append(", Variable name: " + variable_name));
      error_message.append(
        ", Invalid dependency name: " + dependency +
        "\n(HINT: Dependencies that contain \"change(variable_A)\" can only be "
        "used for the field \"variable_A\")\n");

      Assert(reg_possible_deps.find(dependency) != reg_possible_deps.end() ||
               change_possible_deps.at(index).find(dependency) !=
                 change_possible_deps.at(index).end(),
             dealii::ExcMessage(error_message));
    }
}

void
variableAttributeLoader::variableAttributeLoader::validate_postprocess_variable(
  const std::string           &name,
  const std::set<std::string> &name_list,
  const std::set<std::string> &reg_possible_deps,
  const variableAttributes    &pp_variable,
  unsigned int                 index)
{
  Assert(name_list.find(name) == name_list.end(),
         dealii::ExcMessage("PRISMS-PF Error: Postprocess variable names must be "
                            "different from solution variable names.\nProblem index: " +
                            std::to_string(index) + ", Variable name: " + name + "\n"));
  for (const auto &PP_dep : pp_variable.dependencies_PP)
    {
      std::string error_message =
        "PRISMS-PF Error: Invalid postprocessing RHS dependency.\nProblem index: ";
      error_message.append(std::to_string(index).append(", Variable name: " + name));
      error_message.append(", Invalid dependency name: " + PP_dep + "\n");

      Assert(reg_possible_deps.find(PP_dep) != reg_possible_deps.end(),
             dealii::ExcMessage(error_message));
    }
  Assert(pp_variable.dependencies_LHS.empty(),
         dealii::ExcMessage("PRISMS-PF Error: Postprocess variables can only have RHS "
                            "dependencies.\nProblem index: " +
                            std::to_string(index) + ", Variable name: " + name + "\n"));
  Assert(!pp_variable.nucleating_variable && !pp_variable.need_value_nucleation,
         dealii::ExcMessage("PRISMS-PF Error: Postprocess variables cannot be used for "
                            "nucleation.\nProblem index: " +
                            std::to_string(index) + ", Variable name: " + name + "\n"));
}

void
variableAttributeLoader::validate_attributes()
{
  // Make sure main attributes arent set in pp_attributes
  // Make sure dependencies are all variable names and if there are "change()" deps, they
  // are correctly placed
  const std::set<std::pair<std::string, std::string>> reg_delimiters = {
    {"",      "" },
    {"grad(", ")"},
    {"hess(", ")"}
  };
  const std::set<std::string>           forbidden_names = {"grad(", "hess(", "change("};
  std::set<std::string>                 name_list;
  std::set<std::string>                 pp_name_list;
  std::set<std::string>                 reg_possible_deps;
  std::map<uint, std::set<std::string>> change_possible_deps;

  // Populate the expected variable dependencies and check that variable names are mostly
  // well-formed.
  for (const auto &[index, variable] : attributes)
    {
      name_list.insert(variable.name);

      populate_dependencies(reg_delimiters,
                            variable.name,
                            index,
                            reg_possible_deps,
                            change_possible_deps);

      validate_variable_name(variable.name, forbidden_names, "", index);
    }

  // Validate the postprocessed variables
  for (const auto &[pp_index, pp_variable] : pp_attributes)
    {
      pp_name_list.insert(pp_variable.name);

      validate_variable_name(pp_variable.name,
                             forbidden_names,
                             "(postprocess)",
                             pp_index);

      validate_postprocess_variable(pp_variable.name,
                                    name_list,
                                    reg_possible_deps,
                                    pp_variable,
                                    pp_index);
    }
  // Check dependencies
  for (const auto &[index, variable] : attributes)
    {
      validate_dependencies(variable.dependencies_RHS,
                            "RHS",
                            index,
                            variable.name,
                            reg_possible_deps,
                            change_possible_deps);

      validate_dependencies(variable.dependencies_LHS,
                            "LHS",
                            index,
                            variable.name,
                            reg_possible_deps,
                            change_possible_deps);
    }
}

std::string
variableAttributeLoader::strip_whitespace(const std::string &_text)
{
  std::string text = _text;
  text.erase(std::remove(text.begin(), text.end(), ' '), text.end());
  return text;
}