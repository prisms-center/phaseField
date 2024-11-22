#ifndef VARIABLEATTRIBUTELOADER_H
#define VARIABLEATTRIBUTELOADER_H

#include "varTypeEnums.h"
#include "variableAttributes.h"

#include <map>
#include <string>

using EvalFlags = dealii::EvaluationFlags::EvaluationFlags;

/**
 * \brief Class to manage the variable attributes that the user specifies.
 */
class variableAttributeLoader
{
public:
  /**
   * \brief Constructor. Executes the user-facing functions and constructs the variable
   * attributes.
   */
  variableAttributeLoader();

  /**
   * \brief User-facing method where the variable attributes are set.
   */
  void
  loadVariableAttributes();

  /**
   * \brief User-facing method where the postprocessing variable attributes are set.
   */
  void
  loadPostProcessorVariableAttributes();

  /**
   * \brief Set the name of the variable at `index` to `name`.
   *
   * \param index Index of variable
   * \param name Name of variable at `index`
   */
  void
  set_variable_name(const unsigned int &index, const std::string &name) const;

  /**
   * \brief Set the field type of the variable at `index` to `var_type` where `var_type`
   * can be `SCALAR` or `VECTOR`.
   *
   * \param index Index of variable
   * \param var_type Field type of variable at `index` (`SCALAR` or `VECTOR`).
   */
  void
  set_variable_type(const unsigned int &index, const fieldType &var_type) const;

  /**
   * \brief Set the PDE type of the variable at `index` to `var_eq_type` where
   *`var_eq_type`can be `EXPLICIT_TIME_DEPENDENT`, `IMPLICIT_TIME_DEPENDENT`,
   *`TIME_INDEPENDENT`, `AUXILIARY`.
   *
   * \param index Index of variable
   * \param var_eq_type PDE type of variable at `index`.
   */
  void
  set_variable_equation_type(const unsigned int &index, const PDEType &var_eq_type) const;

  /**
   * \brief Set the dependencies for the value term of the RHS equation of the variable at
   * `index`.
   *
   * \param index Index of variable
   * \param dependencies String containing comma-separated list of dependencies for
   * variable at `index` Hint: "variable, grad(variable), hess(variable)"
   */
  void
  set_dependencies_value_term_RHS(const unsigned int &index,
                                  const std::string  &dependencies);

  /**
   * \brief Set the dependencies for the gradient term of the RHS equation of the variable
   * at `index`.
   *
   * \param index Index of variable
   * \param dependencies String containing comma-separated list of dependencies for
   * variable at `index` Hint: "variable, grad(variable), hess(variable)"
   */
  void
  set_dependencies_gradient_term_RHS(const unsigned int &index,
                                     const std::string  &dependencies);

  /**
   * \brief Set the dependencies for the value term of the LHS equation of the variable at
   * `index`.
   *
   * \param index Index of variable
   * \param dependencies String containing comma-separated list of dependencies for
   * variable at `index` Hint: "variable, grad(variable), hess(variable)"
   */
  void
  set_dependencies_value_term_LHS(const unsigned int &index,
                                  const std::string  &dependencies) const;

  /**
   * \brief Set the dependencies for the gradient term of the LHS equation of the variable
   * at `index`.
   *
   * \param index Index of variable
   * \param dependencies String containing comma-separated list of dependencies for
   * variable at `index` Hint: "variable, grad(variable), hess(variable)"
   */
  void
  set_dependencies_gradient_term_LHS(const unsigned int &index,
                                     const std::string  &dependencies) const;

  /**
   * \brief Flag whether the variable at `index` is needed to calculate the nucleation
   * probability.
   *
   * \param index Index of variable
   * \param flag true: variable is needed, false: variable is not needed.
   */
  void
  set_need_value_nucleation(const unsigned int &index, const bool &flag) const;

  /**
   * \brief Flag whether the variable at `index` is can have a nucleation event.
   *
   * \param index Index of variable
   * \param flag true: variable can nucleate, false: variable can not nucleate.
   */
  void
  set_allowed_to_nucleate(const unsigned int &index, const bool &flag) const;

  /**
   * \brief (Postprocess only) Flag whether the postprocessing variable at `index` should
   * have its domain integral calculated and output.
   *
   * \param index Index of variable
   * \param flag true: calculate and output the integral of the field over the domain,
   * false: do nothing
   */
  void
  set_output_integral(const unsigned int &index, const bool &flag) const;

  /**
   * \brief The solutions variable attributes
   */
  std::map<uint, variableAttributes> attributes;

  /**
   * \brief The postprocessing variable attributes
   */
  std::map<uint, variableAttributes> pp_attributes;

  /**
   * \brief Useful pointer for setting whether solution or postprocessiong variables are
   * being loaded
   */
  std::map<uint, variableAttributes> *relevant_attributes = nullptr;

private:
  /**
   * \brief Perform a suite of assertions on attributes and pp_attributes to ensure that
   * the user's inputs are well-formed
   */
  void
  validate_attributes();

  /**
   * \brief Validate that the variable name is not empty and does not contain any
   * forbidden substrings (names).
   */
  static void
  validate_variable_name(const std::string           &name,
                         const std::set<std::string> &forbidden_names,
                         const std::string           &context,
                         unsigned int                 index);

  /**
   * \brief Populate dependencies that we should expect from the user.
   */
  static void
  populate_dependencies(
    const std::set<std::pair<std::string, std::string>> &reg_delimiters,
    const std::string                                   &variable_name,
    unsigned int                                         index,
    std::set<std::string>                               &reg_possible_deps,
    std::map<uint, std::set<std::string>>               &change_possible_deps);

  /**
   * \brief Validate the dependencies (RHS or LHS) that the user has provided.
   */
  static void
  validate_dependencies(
    const std::set<std::string>                 &dependencies,
    const std::string                           &context,
    unsigned int                                 index,
    const std::string                           &variable_name,
    const std::set<std::string>                 &reg_possible_deps,
    const std::map<uint, std::set<std::string>> &change_possible_deps);

  /**
   * \brief Validate the postprocess variables.
   */
  static void
  validate_postprocess_variable(const std::string           &name,
                                const std::set<std::string> &name_list,
                                const std::set<std::string> &reg_possible_deps,
                                const variableAttributes    &pp_variable,
                                unsigned int                 index);

  /**
   * \brief Utility to remove whitespace from strings
   */
  static std::string
  strip_whitespace(const std::string &text);
  // The above function should be moved to a 'utilities' module
};

#endif
