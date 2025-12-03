// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

void
VariableAttributeLoader::init_variable_attributes()
{
  // Load the variable attributes from the user
  load_variable_attributes();

  // Determine the max number of fields that user has defined. This is used to determine
  // the length of the vector for the eval flag set at runtime.
  const auto max_fields = static_cast<Types::Index>(var_attributes.size());

  // Determine the max numer of dependency types. This is used to determine the length of
  // the vector for the eval flag set at runtime.
  const auto max_dependency_types = static_cast<Types::Index>(DependencyType::OldFour);

  // Format the dependencies and add the max fields and dependency types
  for (auto &[index, variable] : var_attributes)
    {
      variable.format_dependencies();
      variable.max_fields           = max_fields;
      variable.max_dependency_types = max_dependency_types;
    }

  // Validate the attributes
  validate_attributes();

  // Parse the string dependencies into the eval flag set
  for (auto &[index, variable] : var_attributes)
    {
      variable.parse_residual_dependencies();
      variable.parse_dependencies(var_attributes, max_fields, max_dependency_types);
    }

  // Validate the old solution dependencies
  validate_old_solution_dependencies();

  // Determine the field solve type
  for (auto &[index, variable] : var_attributes)
    {
      variable.determine_field_solve_type(var_attributes);
    }

  // Compute the set of solve blocks
  std::set<Types::Index> solve_block_collection;
  for (auto &[index, variable] : var_attributes)
    {
      solve_block_collection.insert(variable.solve_block);
    }
  // Compute the shared dependencies
  for (const FieldSolveType field_solve_type : {FieldSolveType::ExplicitConstant,
                                                FieldSolveType::Explicit,
                                                FieldSolveType::ExplicitPostprocess})
    {
      for (const Types::Index solve_block : solve_block_collection)
        {
          auto subset_attributes =
            compute_subset_attributes(var_attributes, field_solve_type, solve_block);
          compute_shared_dependencies(subset_attributes);
        }
    }

  // Print variable attributes to summary.log
  for (const auto &[index, variable] : var_attributes)
    {
      variable.print();
    }
}

std::map<unsigned int, VariableAttributes>
VariableAttributeLoader::get_var_attributes() const
{
  return var_attributes;
}

void
VariableAttributeLoader::set_variable_name(const unsigned int &index,
                                           const std::string  &name)
{
  var_attributes[index].field_index = index;
  var_attributes[index].name        = name;
}

void
VariableAttributeLoader::set_variable_type(const unsigned int          &index,
                                           const FieldInfo::TensorRank &field_type)
{
  switch (field_type)
    {
      case FieldInfo::TensorRank::Scalar:
      case FieldInfo::TensorRank::Vector:
        var_attributes[index].field_info.tensor_rank = field_type;
        break;
      default:
        throw std::invalid_argument(
          "Invalid FieldInfo::TensorRank value provided in set_variable_type().");
    }
}

void
VariableAttributeLoader::set_variable_equation_type(const unsigned int &index,
                                                    const PDEType      &pde_type)
{
  switch (pde_type)
    {
      case PDEType::ExplicitTimeDependent:
      case PDEType::ImplicitTimeDependent:
      case PDEType::TimeIndependent:
      case PDEType::Auxiliary:
      case PDEType::Constant:
        var_attributes[index].pde_type = pde_type;
        break;
      default:
        throw std::invalid_argument(
          "Invalid PDEType value provided in set_variable_equation_type().");
    }
}

void
VariableAttributeLoader::set_is_postprocessed_field(const unsigned int &index,
                                                    const bool         &is_postprocess)
{
  var_attributes[index].is_postprocessed_variable = is_postprocess;
}

template <typename Iterable>
void
VariableAttributeLoader::set_is_nucleation_rate(const unsigned int &index,
                                                const bool         &is_nucleation,
                                                const Iterable     &nucleating_fields)
{
  var_attributes[index].is_nucleation_rate_variable = is_nucleation;
  var_attributes[index].raw_dependencies.nucleating_fields.insert(
    nucleating_fields.begin(),
    nucleating_fields.end());
}

void
VariableAttributeLoader::set_solve_block(const unsigned int &index,
                                         const Types::Index &solve_block)
{
  var_attributes[index].solve_block = solve_block;
}

void
VariableAttributeLoader::set_dependencies_value_term_rhs(const unsigned int &index,
                                                         const std::string  &dependencies)
{
  const std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  insert_dependencies_value_term_rhs(index, dependencies_set);
}

void
VariableAttributeLoader::set_dependencies_gradient_term_rhs(
  const unsigned int &index,
  const std::string  &dependencies)
{
  const std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  insert_dependencies_gradient_term_rhs(index, dependencies_set);
}

void
VariableAttributeLoader::set_dependencies_value_term_lhs(const unsigned int &index,
                                                         const std::string  &dependencies)
{
  const std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  insert_dependencies_value_term_lhs(index, dependencies_set);
}

void
VariableAttributeLoader::set_dependencies_gradient_term_lhs(
  const unsigned int &index,
  const std::string  &dependencies)
{
  const std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  insert_dependencies_gradient_term_lhs(index, dependencies_set);
}

template <typename Iterable>
void
VariableAttributeLoader::insert_dependencies_value_term_rhs(const unsigned int &index,
                                                            const Iterable &dependencies)
{
  var_attributes[index].raw_dependencies.dependencies_value_rhs.insert(
    dependencies.begin(),
    dependencies.end());
}

template <typename Iterable>
void
VariableAttributeLoader::insert_dependencies_gradient_term_rhs(
  const unsigned int &index,
  const Iterable     &dependencies)
{
  var_attributes[index].raw_dependencies.dependencies_gradient_rhs.insert(
    dependencies.begin(),
    dependencies.end());
}

template <typename Iterable>
void
VariableAttributeLoader::insert_dependencies_value_term_lhs(const unsigned int &index,
                                                            const Iterable &dependencies)
{
  var_attributes[index].raw_dependencies.dependencies_value_lhs.insert(
    dependencies.begin(),
    dependencies.end());
}

template <typename Iterable>
void
VariableAttributeLoader::insert_dependencies_gradient_term_lhs(
  const unsigned int &index,
  const Iterable     &dependencies)
{
  var_attributes[index].raw_dependencies.dependencies_gradient_lhs.insert(
    dependencies.begin(),
    dependencies.end());
}

void
VariableAttributeLoader::validate_variable_name(
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
VariableAttributeLoader::populate_dependencies(
  const std::set<std::pair<std::string, std::string>> &reg_delimiters,
  const std::set<std::pair<std::string, std::string>> &dep_type_delimiters,
  const std::string                                   &variable_name,
  unsigned int                                         index,
  std::set<std::string>                               &reg_possible_deps,
  std::map<unsigned int, std::set<std::string>>       &change_possible_deps)
{
  // If we are dealing with a postprocessed variable, no dependencies are valid
  if (var_attributes.at(index).is_postprocessed_variable)
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
VariableAttributeLoader::validate_dependencies(
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

      AssertThrow(reg_possible_deps.contains(dependency) ||
                    change_possible_deps.at(index).contains(dependency),
                  dealii::ExcMessage(error_message));
    }
}

void
VariableAttributeLoader::validate_attributes()
{
  for (const auto &[index, variable] : var_attributes)
    {
      // Check that postprocessed variables are only explicit
      AssertThrow(!(variable.is_postprocessed_variable &&
                    variable.pde_type != PDEType::ExplicitTimeDependent),
                  dealii::ExcMessage(
                    "Currently, postprocessing only allows explicit equations."));
      // Check that nucleation rates are scalars
      AssertThrow(!(variable.is_nucleation_rate_variable &&
                    variable.field_info.tensor_rank != FieldInfo::TensorRank::Scalar),
                  dealii::ExcMessage(
                    "Currently, nucleation rates must be scalar fields."));
      // Check that constant fields have no dependencies
      AssertThrow(!(variable.pde_type == PDEType::Constant) ||
                    (variable.raw_dependencies.dependencies_rhs.empty() &&
                     variable.raw_dependencies.dependencies_lhs.empty()),
                  dealii::ExcMessage("Constant fields are determined by the initial "
                                     "condition. They cannot have dependencies."));
    }

  // Make sure dependencies are all variable names. If there are change() dependencies
  // make sure they are of the same field. Note that this does not enforce them to be on
  // the LHS. It also does not enforce the fact div(), symgrad(), or curl() dependencies
  // must belong to a vector field. Both of these are done later in VariableAttributes.
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
      validate_dependencies(variable.raw_dependencies.dependencies_rhs,
                            "RHS",
                            index,
                            variable.name,
                            reg_possible_deps,
                            change_possible_deps);

      validate_dependencies(variable.raw_dependencies.dependencies_lhs,
                            "LHS",
                            index,
                            variable.name,
                            reg_possible_deps,
                            change_possible_deps);
    }
}

void
VariableAttributeLoader::validate_old_solution_dependencies()
{
  // First create a combined dependency set for all fields
  std::map<std::pair<unsigned int, DependencyType>,
           dealii::EvaluationFlags::EvaluationFlags>
    dependency_set;
  for (const auto &[index, variable] : var_attributes)
    {
      Types::Index field_index = 0;
      for (const auto &local_dependency_set : variable.eval_flag_set_rhs)
        {
          Types::Index dep_index = 0;
          for (const auto &value : local_dependency_set)
            {
              // Skip where the evaluation flags are nothing
              if (value == dealii::EvaluationFlags::EvaluationFlags::nothing)
                {
                  dep_index++;
                  continue;
                }
              dependency_set[std::make_pair(field_index,
                                            static_cast<DependencyType>(dep_index))] |=
                value;
              dep_index++;
            }
          field_index++;
        }

      field_index = 0;
      for (const auto &local_dependency_set : variable.eval_flag_set_lhs)
        {
          Types::Index dep_index = 0;
          for (const auto &value : local_dependency_set)
            {
              // Skip where the evaluation flags are nothing
              if (value == dealii::EvaluationFlags::EvaluationFlags::nothing)
                {
                  dep_index++;
                  continue;
                }
              dependency_set[std::make_pair(field_index,
                                            static_cast<DependencyType>(dep_index))] |=
                value;
              dep_index++;
            }
          field_index++;
        }
    }

  // Check that constant equations do not have old dependencies
  for (const auto &[pair, flag] : dependency_set)
    {
      if (var_attributes.at(pair.first).pde_type == PDEType::Constant)
        {
          AssertThrow(
            pair.second == DependencyType::Normal,
            dealii::ExcMessage(
              "Constant fields cannot be specified as old() or change() dependencies"));
        }
    }

  // Check that old fields dependencies are sequential
  for (const auto &[index, variable] : var_attributes)
    {
      std::vector<DependencyType> old_types = {DependencyType::OldOne,
                                               DependencyType::OldTwo,
                                               DependencyType::OldThree,
                                               DependencyType::OldFour};

      // Find first gap in sequence
      unsigned int gap_index = 0;
      while (gap_index < old_types.size() &&
             dependency_set.contains(std::make_pair(gap_index, old_types[gap_index])))
        {
          gap_index++;
        }

      // Check no dependencies exist after gap
      for (unsigned int second_index = gap_index; second_index < old_types.size();
           second_index++)
        {
          if (dependency_set.contains(std::make_pair(index, old_types[second_index])))
            {
              AssertThrow(false,
                          dealii::ExcMessage(
                            "If old_n() of a field is specified, the "
                            "previous old_n() to old_1() must be present."));
            }
        }
    }
}

std::map<Types::Index, VariableAttributes *>
VariableAttributeLoader::compute_subset_attributes(
  std::map<Types::Index, VariableAttributes> &variable_attributes,
  FieldSolveType                              field_solve_type,
  Types::Index                                solve_priority) const
{
  std::map<Types::Index, VariableAttributes *> local_subset_attributes;

  for (auto &[index, variable] : variable_attributes)
    {
      if (variable.field_solve_type == field_solve_type &&
          variable.solve_block == solve_priority)
        {
          local_subset_attributes.emplace(index, &variable);
        }
    }

  return local_subset_attributes;
}

void
VariableAttributeLoader::compute_shared_dependencies(
  std::map<Types::Index, VariableAttributes *> &variable_attributes)
{
  // If the map is entry return early
  if (variable_attributes.empty())
    {
      return;
    }

  // Grab the first entry in the variable_attributes to get some useful information
  const auto *first_variable = variable_attributes.begin()->second;
  [[maybe_unused]] const FieldSolveType field_solve_type =
    first_variable->field_solve_type;
  const Types::Index max_fields       = first_variable->max_fields;
  const Types::Index max_dependencies = first_variable->max_dependency_types;

  // If the field_solve_type is not something we would expect thow an assertion. We could
  // return early, but I slightly favor an assertion and choosing the appropriate
  // FieldSolveType's upstream.
  Assert(field_solve_type == FieldSolveType::Explicit ||
           field_solve_type == FieldSolveType::ExplicitConstant ||
           field_solve_type == FieldSolveType::ExplicitPostprocess,
         dealii::ExcMessage(
           "compute_shared_dependencies() should only be used for concurrent solves."));

  // Create a vector for the shared eval flags and dependencies
  std::vector<std::vector<dealii::EvaluationFlags::EvaluationFlags>> shared_eval_flags(
    max_fields,
    std::vector<dealii::EvaluationFlags::EvaluationFlags>(
      max_dependencies,
      dealii::EvaluationFlags::EvaluationFlags::nothing));
  std::vector<std::vector<FieldInfo::TensorRank>> shared_dependencies(
    max_fields,
    std::vector<FieldInfo::TensorRank>(max_dependencies,
                                       FieldInfo::TensorRank::Undefined));

  // Populate the shared eval flags
  for (const auto &[index, variable] : variable_attributes)
    {
      for (Types::Index field_index = 0; field_index < max_fields; field_index++)
        {
          for (Types::Index dependency_index = 0; dependency_index < max_dependencies;
               dependency_index++)
            {
              const auto &eval_flag =
                variable->eval_flag_set_rhs.at(field_index).at(dependency_index);
              const auto &dependency_type =
                variable->dependency_set_rhs.at(field_index).at(dependency_index);

              // If the eval flag is nothing and the dependency type is invalid skip it
              if (eval_flag == dealii::EvaluationFlags::EvaluationFlags::nothing &&
                  dependency_type == FieldInfo::TensorRank::Undefined)
                {
                  continue;
                }

              // Check that change terms are not found in the RHS of the dependency set
              Assert(eval_flag == dealii::EvaluationFlags::EvaluationFlags::nothing ||
                       dependency_index !=
                         static_cast<Types::Index>(DependencyType::Change),
                     dealii::ExcMessage(
                       "Change terms are not allowed in the RHS of the dependency set"));

              // If the dependency is not a change term, then we can add the dependency
              shared_eval_flags.at(field_index).at(dependency_index) |= eval_flag;

              // Also add the field type to the shared dependencies
              shared_dependencies.at(field_index).at(dependency_index) =
                variable->dependency_set_rhs.at(field_index).at(dependency_index);
            }
        }
    }

  // Assign the shared dependencies to the subset attributes
  for (auto &[index, variable] : variable_attributes)
    {
      variable->eval_flag_set_rhs  = shared_eval_flags;
      variable->dependency_set_rhs = shared_dependencies;
    }
}

#include "core/variable_attribute_loader.inst"

PRISMS_PF_END_NAMESPACE
