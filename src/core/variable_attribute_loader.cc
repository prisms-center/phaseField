// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <boost/range/algorithm/remove.hpp>

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_attributes.h>

#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

void
variableAttributeLoader::loadVariableAttributes()
{}

void
variableAttributeLoader::init_variable_attributes()
{
  loadVariableAttributes();

  for (auto &[index, variable] : var_attributes)
    {
      variable.format_dependencies();
    }
  validate_attributes();
  for (auto &[index, variable] : var_attributes)
    {
      variable.parse_residual_dependencies();
      variable.parse_dependencies(var_attributes);
    }
  validate_old_solution_dependencies();
  for (auto &[index, variable] : var_attributes)
    {
      variable.determine_field_solve_type(var_attributes);
    }

  // Print variable attributes to summary.log
  for (const auto &[index, variable] : var_attributes)
    {
      variable.print();
    }
}

std::map<unsigned int, variableAttributes>
variableAttributeLoader::get_var_attributes() const
{
  return var_attributes;
}

void
variableAttributeLoader::set_variable_name(const unsigned int &index,
                                           const std::string  &name)
{
  var_attributes[index].field_index = index;
  var_attributes[index].name        = name;
}

void
variableAttributeLoader::set_variable_type(const unsigned int &index,
                                           const fieldType    &field_type)
{
  switch (field_type)
    {
      case fieldType::SCALAR:
      case fieldType::VECTOR:
        var_attributes[index].field_type = field_type;
        break;
      default:
        throw std::invalid_argument(
          "Invalid fieldType value provided in set_variable_type().");
    }
}

void
variableAttributeLoader::set_variable_equation_type(const unsigned int &index,
                                                    const PDEType      &pde_type)
{
  switch (pde_type)
    {
      case PDEType::EXPLICIT_TIME_DEPENDENT:
      case PDEType::IMPLICIT_TIME_DEPENDENT:
      case PDEType::TIME_INDEPENDENT:
      case PDEType::AUXILIARY:
      case PDEType::CONSTANT:
        var_attributes[index].pde_type = pde_type;
        break;
      default:
        throw std::invalid_argument(
          "Invalid PDEType value provided in set_variable_equation_type().");
    }
}

void
variableAttributeLoader::set_is_postprocessed_field(const unsigned int &index,
                                                    const bool         &is_postprocess)
{
  var_attributes[index].is_postprocess = is_postprocess;
}

void
variableAttributeLoader::set_dependencies_value_term_RHS(const unsigned int &index,
                                                         const std::string  &dependencies)
{
  const std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  insert_dependencies_value_term_RHS(index, dependencies_set);
}

void
variableAttributeLoader::set_dependencies_gradient_term_RHS(
  const unsigned int &index,
  const std::string  &dependencies)
{
  const std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  insert_dependencies_gradient_term_RHS(index, dependencies_set);
}

void
variableAttributeLoader::set_dependencies_value_term_LHS(const unsigned int &index,
                                                         const std::string  &dependencies)
{
  const std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  insert_dependencies_value_term_LHS(index, dependencies_set);
}

void
variableAttributeLoader::set_dependencies_gradient_term_LHS(
  const unsigned int &index,
  const std::string  &dependencies)
{
  const std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  insert_dependencies_gradient_term_LHS(index, dependencies_set);
}

template <typename Iterable>
void
variableAttributeLoader::insert_dependencies_value_term_RHS(const unsigned int &index,
                                                            const Iterable &dependencies)
{
  var_attributes[index].dependencies_value_RHS.insert(dependencies.begin(),
                                                      dependencies.end());
}

template <typename Iterable>
void
variableAttributeLoader::insert_dependencies_gradient_term_RHS(
  const unsigned int &index,
  const Iterable     &dependencies)
{
  var_attributes[index].dependencies_gradient_RHS.insert(dependencies.begin(),
                                                         dependencies.end());
}

template <typename Iterable>
void
variableAttributeLoader::insert_dependencies_value_term_LHS(const unsigned int &index,
                                                            const Iterable &dependencies)
{
  var_attributes[index].dependencies_value_LHS.insert(dependencies.begin(),
                                                      dependencies.end());
}

template <typename Iterable>
void
variableAttributeLoader::insert_dependencies_gradient_term_LHS(
  const unsigned int &index,
  const Iterable     &dependencies)
{
  var_attributes[index].dependencies_gradient_LHS.insert(dependencies.begin(),
                                                         dependencies.end());
}

void
variableAttributeLoader::validate_variable_name(
  const std::string           &name,
  const std::set<std::string> &forbidden_names,
  const std::string           &context,
  unsigned int                 index)
{
  AssertThrow(!name.empty(),
              dealii::ExcMessage(context +
                                 " Variable names must not be empty.\nProblem index: " +
                                 std::to_string(index)));

  for ([[maybe_unused]] const std::string &forbidden_name : forbidden_names)
    {
      std::string error_message = context +
                                  " Variable names must not contain \"grad(\", "
                                  "\"hess(\", \"change(\", \"hessdiag(\", \"lap(\", "
                                  "\"div(\", \"symgrad(\", \"curl(\", \"old_1(\", "
                                  "\"old_2(\", \"old_3(\", \"old_4(\".\nProblem index: ";
      error_message.append(std::to_string(index).append(", Variable name: " + name));

      AssertThrow(name.find(forbidden_name) == std::string::npos,
                  dealii::ExcMessage(error_message));
    }
}

void
variableAttributeLoader::populate_dependencies(
  const std::set<std::pair<std::string, std::string>> &reg_delimiters,
  const std::set<std::pair<std::string, std::string>> &dep_type_delimiters,
  const std::string                                   &variable_name,
  unsigned int                                         index,
  std::set<std::string>                               &reg_possible_deps,
  std::map<unsigned int, std::set<std::string>>       &change_possible_deps)
{
  // If we are dealing with a postprocessed variable, no dependencies are valid
  if (var_attributes.at(index).is_postprocess)
    {
      return;
    }

  // Populate dependencies for main variables
  for (const auto &reg_delims : reg_delimiters)
    {
      for (const auto &type_delims : dep_type_delimiters)
        {
          if (type_delims.first == "change(")
            {
              change_possible_deps[index].insert(reg_delims.first + type_delims.first +
                                                 variable_name + type_delims.second +
                                                 reg_delims.second);
            }
          else
            {
              reg_possible_deps.insert(reg_delims.first + type_delims.first +
                                       variable_name + type_delims.second +
                                       reg_delims.second);
            }
        }
    }
}

void
variableAttributeLoader::validate_dependencies(
  const std::set<std::string>                  &dependencies,
  const std::string                            &context,
  unsigned int                                  index,
  const std::string                            &variable_name,
  [[maybe_unused]] const std::set<std::string> &reg_possible_deps,
  [[maybe_unused]] const std::map<unsigned int, std::set<std::string>>
    &change_possible_deps)
{
  for (const std::string &dependency : dependencies)
    {
      std::string error_message = "Invalid " + context + " dependency.\nProblem index: ";
      error_message.append(
        std::to_string(index).append(", Variable name: " + variable_name));
      error_message.append(
        ", Invalid dependency name: " + dependency +
        "\n(HINT: Dependencies that contain \"change(variable_A)\" can only be "
        "used for the field \"variable_A\"). Additionally, postprocessed fields are not "
        "allowed to be dependencies.\n");

      AssertThrow(reg_possible_deps.find(dependency) != reg_possible_deps.end() ||
                    change_possible_deps.at(index).find(dependency) !=
                      change_possible_deps.at(index).end(),
                  dealii::ExcMessage(error_message));
    }
}

void
variableAttributeLoader::validate_attributes()
{
  for (const auto &[index, variable] : var_attributes)
    {
      // Check that postprocessed variables are only explicit
      AssertThrow(
        !variable.is_postprocess || variable.pde_type == PDEType::EXPLICIT_TIME_DEPENDENT,
        dealii::ExcMessage("Currently, postprocessing only allows explicit equations."));
      // Check that constant fields have no dependencies
      AssertThrow(!(variable.pde_type == PDEType::CONSTANT) ||
                    (variable.dependencies_RHS.empty() &&
                     variable.dependencies_LHS.empty()),
                  dealii::ExcMessage("Constant fields are determined by the initial "
                                     "condition. They cannot have dependencies."));
    }

  // Make sure dependencies are all variable names. If there are change() dependencies
  // make sure they are of the same field. Note that this does not enforce them to be on
  // the LHS. It also does not enforce the fact div(), symgrad(), or curl() dependencies
  // must belong to a vector field. Both of these are done later in variableAttributes.
  const std::set<std::pair<std::string, std::string>> reg_delimiters = {
    {"",          "" },
    {"grad(",     ")"},
    {"hess(",     ")"},
    {"hessdiag(", ")"},
    {"lap(",      ")"},
    {"div(",      ")"},
    {"symgrad(",  ")"},
    {"curl(",     ")"},
  };

  const std::set<std::pair<std::string, std::string>> dep_type_delimiters = {
    {"",        "" },
    {"change(", ")"},
    {"old_1(",  ")"},
    {"old_2(",  ")"},
    {"old_3(",  ")"},
    {"old_4(",  ")"},
  };

  const std::set<std::string> forbidden_names = {"grad(",
                                                 "hess(",
                                                 "change(",
                                                 "hessdiag(",
                                                 "lap(",
                                                 "div(",
                                                 "symgrad(",
                                                 "curl(",
                                                 "old_1(",
                                                 "old_2(",
                                                 "old_3(",
                                                 "old_4("};

  std::set<std::string>                         name_list;
  std::set<std::string>                         reg_possible_deps;
  std::map<unsigned int, std::set<std::string>> change_possible_deps;

  // Populate the expected variable dependencies and check that variable names are mostly
  // well-formed.
  for (const auto &[index, variable] : var_attributes)
    {
      name_list.insert(variable.name);

      populate_dependencies(reg_delimiters,
                            dep_type_delimiters,
                            variable.name,
                            index,
                            reg_possible_deps,
                            change_possible_deps);

      validate_variable_name(variable.name, forbidden_names, "", index);
    }

  // Check dependencies
  for (const auto &[index, variable] : var_attributes)
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

void
variableAttributeLoader::validate_old_solution_dependencies()
{
  // First create a combined dependency set
  std::unordered_map<std::pair<unsigned int, dependencyType>,
                     dealii::EvaluationFlags::EvaluationFlags,
                     pairHash>
    dependency_set;
  for (const auto &[index, variable] : var_attributes)
    {
      if (!variable.eval_flag_set_RHS.empty())
        {
          for (const auto &[pair, flag] : variable.eval_flag_set_RHS)
            {
              dependency_set[pair] |= flag;
            }
        }
      if (!variable.eval_flag_set_LHS.empty())
        {
          for (const auto &[pair, flag] : variable.eval_flag_set_LHS)
            {
              dependency_set[pair] |= flag;
            }
        }
    }

  // Check that constant equations do not have old dependencies
  for (const auto &[pair, flag] : dependency_set)
    {
      if (var_attributes.at(pair.first).pde_type == PDEType::CONSTANT)
        {
          AssertThrow(
            pair.second == dependencyType::NORMAL,
            dealii::ExcMessage(
              "Constant fields cannot be specified as old() or change() dependencies"));
        }
    }

  // Check that old fields dependencies are sequential
  for (const auto &[index, variable] : var_attributes)
    {
      const auto old_1 = std::make_pair(index, dependencyType::OLD_1);
      const auto old_2 = std::make_pair(index, dependencyType::OLD_2);
      const auto old_3 = std::make_pair(index, dependencyType::OLD_3);
      const auto old_4 = std::make_pair(index, dependencyType::OLD_4);

      if (dependency_set.find(old_1) != dependency_set.end())
        {
          if (dependency_set.find(old_2) != dependency_set.end())
            {
              if (dependency_set.find(old_3) != dependency_set.end())
                {
                  if (dependency_set.find(old_4) != dependency_set.end())
                    {
                      return;
                    }
                }
              else if (dependency_set.find(old_4) != dependency_set.end())
                {
                  AssertThrow(false,
                              dealii::ExcMessage(
                                "If old_n() of a field is specified, the "
                                "previous old_n() to old_1() must be present."));
                }
              else
                {
                  return;
                }
            }
          else if (dependency_set.find(old_3) != dependency_set.end() ||
                   dependency_set.find(old_4) != dependency_set.end())
            {
              AssertThrow(false,
                          dealii::ExcMessage(
                            "If old_n() of a field is specified, the "
                            "previous old_n() to old_1() must be present."));
            }
          else
            {
              return;
            }
        }
      else if (dependency_set.find(old_2) != dependency_set.end() ||
               dependency_set.find(old_3) != dependency_set.end() ||
               dependency_set.find(old_4) != dependency_set.end())
        {
          AssertThrow(false,
                      dealii::ExcMessage("If old_n() of a field is specified, the "
                                         "previous old_n() to old_1() must be present."));
        }
      else
        {
          return;
        }

      AssertThrow(false,
                  dealii::ExcMessage("If old_n() of a field is specified, the "
                                     "previous old_n() to old_1() must be present."));
    }
}

std::string
variableAttributeLoader::strip_whitespace(const std::string &_text)
{
  std::string text = _text;
  text.erase(boost::range::remove(text, ' '), text.end());
  return text;
}

PRISMS_PF_END_NAMESPACE
