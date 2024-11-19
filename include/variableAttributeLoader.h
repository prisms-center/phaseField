// Class to manage the variable attributes that the user specifies
#ifndef VARIABLEATTRIBUTELOADER_H
#define VARIABLEATTRIBUTELOADER_H

#include "varTypeEnums.h"
#include "variableAttributes.h"

#include <map>
#include <string>

using EvalFlags = dealii::EvaluationFlags::EvaluationFlags;

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
  void
  set_variable_name(const unsigned int &index, const std::string &name);

  void
  set_variable_type(const unsigned int &index, const fieldType &);

  void
  set_variable_equation_type(const unsigned int &index, const PDEType &);

  void
  set_dependencies_value_term_RHS(const unsigned int &index,
                                  const std::string  &dependencies);

  void
  set_dependencies_gradient_term_RHS(const unsigned int &index,
                                     const std::string  &dependencies);

  void
  set_dependencies_value_term_LHS(const unsigned int &index,
                                  const std::string  &dependencies);

  void
  set_dependencies_gradient_term_LHS(const unsigned int &index,
                                     const std::string  &dependencies);

  void
  set_need_value_nucleation(const unsigned int &index, const bool &);

  void
  set_allowed_to_nucleate(const unsigned int &index, const bool &);

  void
  set_output_integral(const unsigned int &index, const bool &);

  void
  format_dependencies();
  void
  validate_attributes();
  std::string
  strip_whitespace(const std::string &text);

  // Members
  std::map<uint, variableAttributes>  attributes;
  std::map<uint, variableAttributes>  pp_attributes;
  std::map<uint, variableAttributes> *relevant_attributes = nullptr;

  unsigned int number_of_variables;
  unsigned int pp_number_of_variables;
};

#endif
