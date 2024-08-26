// Class to hold the variable attributes that will be passed to a
// userInputParameters object
#ifndef VARIABLEATTRIBUTELOADER_H
#define VARIABLEATTRIBUTELOADER_H

#include "EquationDependencyParser.h"
#include "model_variables.h"
#include "varTypeEnums.h"

#include <string>
#include <vector>

class variableAttributeLoader
{
public:
  // Constructor
  variableAttributeLoader();

  // Methods where the attributes are set
  void
  loadVariableAttributes();
  void
  loadPostProcessorVariableAttributes();

  // Methods to set the parameter_attributes
  bool setting_primary_field_attributes;
  void
  set_variable_name(unsigned int index, std::string name);
  void
  set_variable_type(unsigned int index, fieldType);
  void
  set_variable_equation_type(unsigned int index, PDEType);

  void
  set_need_value_nucleation(unsigned int index, bool);
  void
  set_allowed_to_nucleate(unsigned int index, bool);

  void
  set_output_integral(unsigned int index, bool);

  // Variable inputs (v2.0)
  std::vector<std::pair<unsigned int, std::string>> var_name_list;
  std::vector<std::pair<unsigned int, fieldType>>   var_type_list;
  std::vector<std::pair<unsigned int, PDEType>>     var_eq_type_list;
  std::vector<std::pair<unsigned int, bool>>        need_value_list;
  std::vector<std::pair<unsigned int, bool>>        need_gradient_list;
  std::vector<std::pair<unsigned int, bool>>        need_hessian_list;
  std::vector<std::pair<unsigned int, bool>>        need_value_residual_list;
  std::vector<std::pair<unsigned int, bool>>        need_gradient_residual_list;
  std::vector<std::pair<unsigned int, bool>>        need_value_list_LHS;
  std::vector<std::pair<unsigned int, bool>>        need_gradient_list_LHS;
  std::vector<std::pair<unsigned int, bool>>        need_hessian_list_LHS;
  std::vector<std::pair<unsigned int, bool>>        need_value_residual_list_LHS;
  std::vector<std::pair<unsigned int, bool>>        need_gradient_residual_list_LHS;

  std::vector<std::pair<unsigned int, bool>> need_value_change_list_LHS;
  std::vector<std::pair<unsigned int, bool>> need_gradient_change_list_LHS;
  std::vector<std::pair<unsigned int, bool>> need_hessian_change_list_LHS;

  std::vector<std::pair<unsigned int, bool>> need_value_list_PP;
  std::vector<std::pair<unsigned int, bool>> need_gradient_list_PP;
  std::vector<std::pair<unsigned int, bool>> need_hessian_list_PP;
  std::vector<std::pair<unsigned int, bool>> need_value_list_nucleation;
  std::vector<std::pair<unsigned int, bool>> nucleating_variable_list;

  std::vector<std::pair<unsigned int, std::string>> var_name_list_PP;
  std::vector<std::pair<unsigned int, fieldType>>   var_type_list_PP;
  std::vector<std::pair<unsigned int, bool>>        output_integral_list;
  std::vector<std::pair<unsigned int, bool>>        need_value_residual_list_PP;
  std::vector<std::pair<unsigned int, bool>>        need_gradient_residual_list_PP;

  void
  set_dependencies_value_term_RHS(unsigned int index, std::string dependencies);
  void
  set_dependencies_gradient_term_RHS(unsigned int index, std::string dependencies);
  void
  set_dependencies_value_term_LHS(unsigned int index, std::string dependencies);
  void
  set_dependencies_gradient_term_LHS(unsigned int index, std::string dependencies);

  std::vector<std::pair<unsigned int, std::string>> var_eq_dependencies_value_RHS;
  std::vector<std::pair<unsigned int, std::string>> var_eq_dependencies_gradient_RHS;
  std::vector<std::pair<unsigned int, std::string>> var_eq_dependencies_value_LHS;
  std::vector<std::pair<unsigned int, std::string>> var_eq_dependencies_gradient_LHS;

  std::vector<std::pair<unsigned int, std::string>> var_eq_dependencies_value_PP;
  std::vector<std::pair<unsigned int, std::string>> var_eq_dependencies_gradient_PP;

  EquationDependencyParser equation_dependency_parser;

  unsigned int number_of_variables;

  std::vector<std::string> var_name;
  std::vector<fieldType>   var_type;
  std::vector<PDEType>     var_eq_type;

  std::vector<bool> var_nonlinear;

  std::vector<bool> nucleating_variable;
  std::vector<bool> need_value_nucleation;

  unsigned int pp_number_of_variables;

  std::vector<std::string> pp_var_name;
  std::vector<fieldType>   pp_var_type;
  std::vector<bool>        pp_calc_integral;
};

#endif
