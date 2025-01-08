#include <deal.II/base/utilities.h>

#include <boost/range/algorithm/remove.hpp>

#include <core/exceptions.h>
#include <core/variable_attribute_loader.h>

void
variableAttributeLoader::loadVariableAttributes()
{}

void
variableAttributeLoader::loadPostProcessorVariableAttributes()
{}

void
variableAttributeLoader::init_variable_attributes()
{
  relevant_attributes = &var_attributes;
  loadVariableAttributes();
  relevant_attributes = &pp_attributes;
  loadPostProcessorVariableAttributes();
  relevant_attributes = nullptr;

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

AttributesList
variableAttributeLoader::get_var_attributes() const
{
  return var_attributes;
}

AttributesList
variableAttributeLoader::get_pp_attributes() const
{
  return pp_attributes;
}

// Methods to set the various variable attributes
void
variableAttributeLoader::set_variable_name(const uint        &index,
                                           const std::string &name) const
{
  (*relevant_attributes)[index].name = name;

  // Also set the field index when setting the name. This should be placed somewhere else
  // in the future.
  (*relevant_attributes)[index].field_index = index;
}

void
variableAttributeLoader::set_variable_type(const uint      &index,
                                           const fieldType &field_type) const
{
  (*relevant_attributes)[index].field_type = field_type;
}

void
variableAttributeLoader::set_variable_equation_type(const uint    &index,
                                                    const PDEType &pde_type) const
{
  (*relevant_attributes)[index].pde_type = pde_type;
}

void
variableAttributeLoader::set_dependencies_value_term_RHS(const uint        &index,
                                                         const std::string &dependencies)
{
  const std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  insert_dependencies_value_term_RHS(index, dependencies_set);
}

void
variableAttributeLoader::set_dependencies_gradient_term_RHS(
  const uint        &index,
  const std::string &dependencies)
{
  const std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  insert_dependencies_gradient_term_RHS(index, dependencies_set);
}

void
variableAttributeLoader::set_dependencies_value_term_LHS(
  const uint        &index,
  const std::string &dependencies) const
{
  const std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  insert_dependencies_value_term_LHS(index, dependencies_set);
}

void
variableAttributeLoader::set_dependencies_gradient_term_LHS(
  const uint        &index,
  const std::string &dependencies) const
{
  const std::vector<std::string> dependencies_set =
    dealii::Utilities::split_string_list(strip_whitespace(dependencies));
  insert_dependencies_gradient_term_LHS(index, dependencies_set);
}

template <typename Iterable>
void
variableAttributeLoader::insert_dependencies_value_term_RHS(const uint     &index,
                                                            const Iterable &dependencies)
{
  /* (*relevant_attributes)[index].dependencies_value_RHS.insert(dependencies.begin(),
   * dependencies.end()); */
  if (relevant_attributes != &pp_attributes)
    {
      var_attributes[index].dependencies_value_RHS.insert(dependencies.begin(),
                                                          dependencies.end());
    }
}

template <typename Iterable>
void
variableAttributeLoader::insert_dependencies_gradient_term_RHS(
  const uint     &index,
  const Iterable &dependencies)
{
  /* (*relevant_attributes)[index].dependencies_gradient_RHS.insert(dependencies.begin(),
   * dependencies.end()); */
  if (relevant_attributes != &pp_attributes)
    {
      var_attributes[index].dependencies_gradient_RHS.insert(dependencies.begin(),
                                                             dependencies.end());
    }
}

template <typename Iterable>
void
variableAttributeLoader::insert_dependencies_value_term_LHS(
  const uint     &index,
  const Iterable &dependencies) const
{
  (*relevant_attributes)[index].dependencies_value_LHS.insert(dependencies.begin(),
                                                              dependencies.end());
}

template <typename Iterable>
void
variableAttributeLoader::insert_dependencies_gradient_term_LHS(
  const uint     &index,
  const Iterable &dependencies) const
{
  (*relevant_attributes)[index].dependencies_gradient_LHS.insert(dependencies.begin(),
                                                                 dependencies.end());
}

void
variableAttributeLoader::validate_variable_name(
  const std::string           &name,
  const std::set<std::string> &forbidden_names,
  const std::string           &context,
  uint                         index)
{
  Assert(!name.empty(),
         dealii::ExcMessage("PRISMS-PF Error: " + context +
                            " Variable names must not be empty.\nProblem index: " +
                            std::to_string(index)));

  for ([[maybe_unused]] const std::string &forbidden_name : forbidden_names)
    {
      std::string error_message = "PRISMS-PF Error: " + context +
                                  " Variable names must not contain \"grad(\", "
                                  "\"hess(\", \"change(\", \"hessdiag(\", \"lap(\", "
                                  "\"div(\", \"symgrad(\", \"curl(\", \"old_1(\", "
                                  "\"old_2(\", \"old_3(\", \"old_4(\".\nProblem index: ";
      error_message.append(std::to_string(index).append(", Variable name: " + name));

      Assert(name.find(forbidden_name) == std::string::npos,
             dealii::ExcMessage(error_message));
    }
}

void
variableAttributeLoader::populate_dependencies(
  const std::set<std::pair<std::string, std::string>> &reg_delimiters,
  const std::set<std::pair<std::string, std::string>> &dep_type_delimiters,
  const std::string                                   &variable_name,
  uint                                                 index,
  std::set<std::string>                               &reg_possible_deps,
  std::map<uint, std::set<std::string>>               &change_possible_deps)
{
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
  const std::set<std::string>                                  &dependencies,
  const std::string                                            &context,
  uint                                                          index,
  const std::string                                            &variable_name,
  [[maybe_unused]] const std::set<std::string>                 &reg_possible_deps,
  [[maybe_unused]] const std::map<uint, std::set<std::string>> &change_possible_deps)
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
variableAttributeLoader::validate_attributes()
{
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

  std::set<std::string>                 name_list;
  std::set<std::string>                 reg_possible_deps;
  std::map<uint, std::set<std::string>> change_possible_deps;

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

std::string
variableAttributeLoader::strip_whitespace(const std::string &_text)
{
  std::string text = _text;
  text.erase(boost::range::remove(text, ' '), text.end());
  return text;
}