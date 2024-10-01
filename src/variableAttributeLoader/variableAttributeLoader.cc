#include "../../include/variableAttributeLoader.h"

#include "../../include/sortIndexEntryPairList.h"

// Constructor
variableAttributeLoader::variableAttributeLoader()
{
  setting_primary_field_attributes = true;
  loadVariableAttributes(); // This is the user-facing function

  number_of_variables = var_name_list.size();

  var_name    = sortIndexEntryPairList(var_name_list, number_of_variables, "var");
  var_type    = sortIndexEntryPairList(var_type_list, number_of_variables, SCALAR);
  var_eq_type = sortIndexEntryPairList(var_eq_type_list,
                                       number_of_variables,
                                       EXPLICIT_TIME_DEPENDENT);

  std::vector<std::string> sorted_dependencies_value_RHS =
    sortIndexEntryPairList(var_eq_dependencies_value_RHS, number_of_variables, "");

  std::vector<std::string> sorted_dependencies_gradient_RHS =
    sortIndexEntryPairList(var_eq_dependencies_gradient_RHS, number_of_variables, "");

  std::vector<std::string> sorted_dependencies_value_LHS =
    sortIndexEntryPairList(var_eq_dependencies_value_LHS, number_of_variables, "");

  std::vector<std::string> sorted_dependencies_gradient_LHS =
    sortIndexEntryPairList(var_eq_dependencies_gradient_LHS, number_of_variables, "");

  nucleating_variable =
    sortIndexEntryPairList(nucleating_variable_list, number_of_variables, false);
  need_value_nucleation =
    sortIndexEntryPairList(need_value_list_nucleation, number_of_variables, false);

  equation_dependency_parser.parse(var_name,
                                   var_eq_type,
                                   sorted_dependencies_value_RHS,
                                   sorted_dependencies_gradient_RHS,
                                   sorted_dependencies_value_LHS,
                                   sorted_dependencies_gradient_LHS,
                                   var_nonlinear);

  setting_primary_field_attributes = false;
  loadPostProcessorVariableAttributes(); // This is the user-facing function

  pp_number_of_variables = var_name_list_PP.size();

  pp_var_name = sortIndexEntryPairList(var_name_list_PP, pp_number_of_variables, "var");
  pp_var_type = sortIndexEntryPairList(var_type_list_PP, pp_number_of_variables, SCALAR);

  std::vector<std::string> pp_sorted_dependencies_value =
    sortIndexEntryPairList(var_eq_dependencies_value_PP, pp_number_of_variables, "");

  std::vector<std::string> pp_sorted_dependencies_gradient =
    sortIndexEntryPairList(var_eq_dependencies_gradient_PP, pp_number_of_variables, "");

  pp_calc_integral =
    sortIndexEntryPairList(output_integral_list, pp_number_of_variables, false);

  equation_dependency_parser.pp_parse(var_name,
                                      pp_var_name,
                                      pp_sorted_dependencies_value,
                                      pp_sorted_dependencies_gradient);
}

// Methods to set the various variable attributes
void
variableAttributeLoader::set_variable_name(unsigned int index, std::string name)
{
  std::pair<unsigned int, std::string> var_pair;
  var_pair.first  = index;
  var_pair.second = std::move(name);

  if (setting_primary_field_attributes)
    {
      var_name_list.push_back(var_pair);
    }
  else
    {
      var_name_list_PP.push_back(var_pair);
    }
}

void
variableAttributeLoader::set_variable_type(unsigned int index, fieldType var_type)
{
  std::pair<unsigned int, fieldType> var_pair;
  var_pair.first  = index;
  var_pair.second = var_type;

  if (setting_primary_field_attributes)
    {
      var_type_list.push_back(var_pair);
    }
  else
    {
      var_type_list_PP.push_back(var_pair);
    }
}

void
variableAttributeLoader::set_variable_equation_type(unsigned int index,
                                                    PDEType      var_eq_type)
{
  std::pair<unsigned int, PDEType> var_pair;
  var_pair.first  = index;
  var_pair.second = var_eq_type;
  var_eq_type_list.push_back(var_pair);
}

void
variableAttributeLoader::set_need_value_nucleation(unsigned int index, bool flag)
{
  std::pair<unsigned int, bool> var_pair;
  var_pair.first  = index;
  var_pair.second = flag;
  need_value_list_nucleation.push_back(var_pair);
}

void
variableAttributeLoader::set_allowed_to_nucleate(unsigned int index, bool flag)
{
  std::pair<unsigned int, bool> var_pair;
  var_pair.first  = index;
  var_pair.second = flag;
  nucleating_variable_list.push_back(var_pair);
}

void
variableAttributeLoader::set_output_integral(unsigned int index, bool flag)
{
  std::pair<unsigned int, bool> var_pair;
  var_pair.first  = index;
  var_pair.second = flag;
  output_integral_list.push_back(var_pair);
}

void
variableAttributeLoader::set_dependencies_value_term_RHS(unsigned int index,
                                                         std::string  dependencies)
{
  std::pair<unsigned int, std::string> var_pair;
  var_pair.first  = index;
  var_pair.second = std::move(dependencies);

  if (setting_primary_field_attributes)
    {
      var_eq_dependencies_value_RHS.push_back(var_pair);
    }
  else
    {
      var_eq_dependencies_value_PP.push_back(var_pair);
    }
}

void
variableAttributeLoader::set_dependencies_gradient_term_RHS(unsigned int index,
                                                            std::string  dependencies)
{
  std::pair<unsigned int, std::string> var_pair;
  var_pair.first  = index;
  var_pair.second = std::move(dependencies);

  if (setting_primary_field_attributes)
    {
      var_eq_dependencies_gradient_RHS.push_back(var_pair);
    }
  else
    {
      var_eq_dependencies_gradient_PP.push_back(var_pair);
    }
}

void
variableAttributeLoader::set_dependencies_value_term_LHS(unsigned int index,
                                                         std::string  dependencies)
{
  std::pair<unsigned int, std::string> var_pair;
  var_pair.first  = index;
  var_pair.second = std::move(dependencies);
  var_eq_dependencies_value_LHS.push_back(var_pair);
}

void
variableAttributeLoader::set_dependencies_gradient_term_LHS(unsigned int index,
                                                            std::string  dependencies)
{
  std::pair<unsigned int, std::string> var_pair;
  var_pair.first  = index;
  var_pair.second = std::move(dependencies);
  var_eq_dependencies_gradient_LHS.push_back(var_pair);
}
