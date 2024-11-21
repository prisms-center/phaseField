#include "../../include/variableAttributeLoader.h"

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>

variableAttributeLoader::variableAttributeLoader()
{
  relevant_attributes = &attributes;
  loadVariableAttributes(); // This is a user-facing function
  relevant_attributes = &pp_attributes;
  loadPostProcessorVariableAttributes(); // This is a user-facing function
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
  const std::set<std::pair<std::string, std::string>> reg_delimiters = {
    {"",      "" },
    {"grad(", ")"},
    {"hess(", ")"}
  };
  const std::set<std::string>           forbidden_names = {"grad(", "hess(", "change("};
  std::set<std::string>                 name_list, pp_name_list, reg_possible_deps;
  std::map<uint, std::set<std::string>> change_possible_deps;
  // Populate name_list, pp_name_list, reg_possible_deps, change_possible_deps
  // Check that variable names are mostly well-formed
  for (const auto &[index, variable] : attributes)
    {
      name_list.insert(variable.name);
      for (const auto &delims : reg_delimiters)
        {
          reg_possible_deps.insert(delims.first + variable.name + delims.second);
          change_possible_deps[index].insert(delims.first + "change(" + variable.name +
                                             ")" + delims.second);
        }
      Assert(!variable.name.empty(),
             dealii::ExcMessage(
               "PRISMS-PF Error: Variable names must not be empty.\nProblem index: " +
               std::to_string(index) + "\n"));
      for ([[maybe_unused]] const std::string &forbidden_name : forbidden_names)
        {
          Assert(variable.name.find(forbidden_name) == std::string::npos,
                 dealii::ExcMessage(
                   "PRISMS-PF Error: Variable names must not contain \"grad()\", "
                   "\"hess()\", \"change()\".\nProblem index: " +
                   std::to_string(index) + ", Variable name: " + variable.name + "\n"));
        }
    }
  for (const auto &[pp_index, pp_variable] : pp_attributes)
    {
      pp_name_list.insert(pp_variable.name);
      Assert(!pp_variable.name.empty(),
             dealii::ExcMessage(
               "PRISMS-PF Error: (postprocess) Variable names must not be "
               "empty.\nProblem index: " +
               std::to_string(pp_index) + "\n"));
      for ([[maybe_unused]] const std::string &forbidden_name : forbidden_names)
        {
          Assert(pp_variable.name.find(forbidden_name) == std::string::npos,
                 dealii::ExcMessage("PRISMS-PF Error: (postprocess) Variable names must "
                                    "not contain \"grad()\", "
                                    "\"hess()\", \"change()\".\nProblem index: " +
                                    std::to_string(pp_index) +
                                    ", Variable name: " + pp_variable.name + "\n"));
        }
    }
  // Check Dependencies & PP
  for (const auto &[index, variable] : attributes)
    {
      for ([[maybe_unused]] const std::string &RHS_dep : variable.dependencies_RHS)
        {
          Assert(
            reg_possible_deps.find(RHS_dep) != reg_possible_deps.end() ||
              change_possible_deps.at(index).find(RHS_dep) !=
                change_possible_deps.at(index).end(),
            dealii::ExcMessage(
              "PRISMS-PF Error: Invalid RHS dependency.\nProblem index: " +
              std::to_string(index) + ", Variable name: " + variable.name +
              ", Invalid dependency name: " + RHS_dep +
              "\n(HINT: Dependencies that contain \"change(variable_A)\" can only be "
              "used for the "
              "field \"variable_A\")\n"));
        }
      for ([[maybe_unused]] const std::string &LHS_dep : variable.dependencies_LHS)
        {
          Assert(
            reg_possible_deps.find(LHS_dep) != reg_possible_deps.end() ||
              change_possible_deps.at(index).find(LHS_dep) !=
                change_possible_deps.at(index).end(),
            dealii::ExcMessage(
              "PRISMS-PF Error: Invalid LHS dependency.\nProblem index: " +
              std::to_string(index) + ", Variable name: " + variable.name +
              ", Invalid dependency name: " + LHS_dep +
              "\n(HINT: Dependencies that contain \"change(variable_A)\" can only be "
              "used for the "
              "field \"variable_A\")\n"));
        }
    }
  for (const auto &[pp_index, pp_variable] : pp_attributes)
    {
      Assert(name_list.find(pp_variable.name) == name_list.end(),
             dealii::ExcMessage(
               "PRISMS-PF Error: Postprocess variable names must be named differently "
               "than solution variable names.\nProblem postprocessing index: " +
               std::to_string(pp_index) + ", Variable name: " + pp_variable.name + "\n"));
      for ([[maybe_unused]] const std::string &PP_dep : pp_variable.dependencies_PP)
        {
          Assert(reg_possible_deps.find(PP_dep) != reg_possible_deps.end(),
                 dealii::ExcMessage(
                   "PRISMS-PF Error: Invalid postprocessing RHS dependency.\nProblem "
                   "postprocessing index: " +
                   std::to_string(pp_index) + ", Variable name: " + pp_variable.name +
                   ", Invalid dependency name: " + PP_dep + "\n"));
        }
      Assert(pp_variable.dependencies_LHS.empty(),
             dealii::ExcMessage(
               "PRISMS-PF Error: Warning: Postprocess variables can only have RHS "
               "dependencies.\nProblem postprocessing index: " +
               std::to_string(pp_index) + ", Variable name: " + pp_variable.name + "\n"));
      Assert(!pp_variable.nucleating_variable && !pp_variable.need_value_nucleation,
             dealii::ExcMessage(
               "PRISMS-PF Error: Warning: Postprocess variables cannot "
               "be used for nucleation.\nProblem postprocessing index: " +
               std::to_string(pp_index) + ", Variable name: " + pp_variable.name + "\n"));
    }
}

std::string
variableAttributeLoader::strip_whitespace(const std::string &_text)
{
  std::string text = _text;
  text.erase(std::remove(text.begin(), text.end(), ' '), text.end());
  return text;
}