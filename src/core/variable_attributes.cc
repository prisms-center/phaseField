#include <boost/range/algorithm/find.hpp>

#include <core/conditional_ostreams.h>
#include <core/variable_attributes.h>
#include <utilities/utilities.h>

void
variableAttributes::format_dependencies()
{
  dependencies_RHS.insert(dependencies_value_RHS.begin(), dependencies_value_RHS.end());
  dependencies_RHS.insert(dependencies_gradient_RHS.begin(),
                          dependencies_gradient_RHS.end());

  dependencies_LHS.insert(dependencies_value_LHS.begin(), dependencies_value_LHS.end());
  dependencies_LHS.insert(dependencies_gradient_LHS.begin(),
                          dependencies_gradient_LHS.end());
}

void
variableAttributes::parse_residual_dependencies()
{
  // Check if either is empty and set value and gradient flags for the
  // residual appropriately. Note that this relies on the fact that only valid
  // dependencies are part of the set.
  if (!dependencies_value_RHS.empty())
    {
      eval_flags_residual_RHS |= dealii::EvaluationFlags::values;
    }
  if (!dependencies_gradient_RHS.empty())
    {
      eval_flags_residual_RHS |= dealii::EvaluationFlags::gradients;
    }
  if (!dependencies_value_LHS.empty())
    {
      eval_flags_residual_LHS |= dealii::EvaluationFlags::values;
    }
  if (!dependencies_gradient_LHS.empty())
    {
      eval_flags_residual_LHS |= dealii::EvaluationFlags::gradients;
    }
}

void
variableAttributes::parse_dependencies(AttributesList &other_var_attributes)
{
  const std::map<std::string, std::pair<dependencyType, EvalFlags>> relevant_flag = []()
  {
    // Modifiers for the dependency types
    std::vector<std::pair<std::string, dependencyType>> dependency_types = {
      {"",        dependencyType::NORMAL},
      {"_change", dependencyType::CHANGE},
      {"_old_1",  dependencyType::OLD_1 },
      {"_old_2",  dependencyType::OLD_2 },
      {"_old_3",  dependencyType::OLD_3 },
      {"_old_4",  dependencyType::OLD_4 },
    };

    // Dependency evaluation types & their corresponding evaluation flag
    std::vector<std::pair<std::string, EvalFlags>> eval_flags = {
      {"value",              dealii::EvaluationFlags::values   },
      {"gradient",           dealii::EvaluationFlags::gradients},
      {"hessian",            dealii::EvaluationFlags::hessians },
      {"hessian_diagonal",   dealii::EvaluationFlags::hessians },
      {"laplacian",          dealii::EvaluationFlags::hessians },
      {"divergence",         dealii::EvaluationFlags::gradients},
      {"symmetric_gradient", dealii::EvaluationFlags::gradients},
      {"curl",               dealii::EvaluationFlags::gradients},
    };

    std::map<std::string, std::pair<dependencyType, EvalFlags>> map;
    for (const auto &dep : dependency_types)
      {
        for (const auto &eval : eval_flags)
          {
            map[eval.first + dep.first] = {dep.second, eval.second};
          }
      }

    return map;
  }();

  const std::map<std::string, std::pair<std::string, std::string>> delimiters = []()
  {
    // Delimiters for the dependency types. Note that there are
    // various overloads for the delimiters.
    std::vector<std::pair<std::string, std::pair<std::string, std::string>>>
      dependency_types_delimiters = {
        {"",        {"", ""}        },
        {"_change", {"change(", ")"}},
        {"_old_1",  {"old_1(", ")"} },
        {"_old_2",  {"old_2(", ")"} },
        {"_old_3",  {"old_3(", ")"} },
        {"_old_4",  {"old_4(", ")"} },
    };

    // Dependency evaluation types & their corresponding delimiter.
    std::vector<std::pair<std::string, std::pair<std::string, std::string>>>
      base_delimiters = {
        {"value",              {"", ""}          },
        {"gradient",           {"grad(", ")"}    },
        {"hessian",            {"hess(", ")"}    },
        {"hessian_diagonal",   {"hessdiag(", ")"}},
        {"laplacian",          {"lap(", ")"}     },
        {"divergence",         {"div(", ")"}     },
        {"symmetric_gradient", {"symgrad(", ")"} },
        {"curl",               {"curl(", ")"}    },
    };

    std::map<std::string, std::pair<std::string, std::string>> map;
    for (const auto &dep : dependency_types_delimiters)
      {
        for (const auto &base : base_delimiters)
          {
            map[base.first + dep.first] = {base.second.first + dep.second.first,
                                           dep.second.second + base.second.second};
          }
      }

    return map;
  }();

  // Helper lambda to validate and fill in dependency sets
  auto set_dependencies =
    [&](const std::set<std::string> &dependencies,
        std::unordered_map<std::pair<uint, dependencyType>, EvalFlags, pairHash>
                          &eval_flag_set,
        const std::string &context)
  {
    // Loop through the available delimiters
    for (const auto &[variation, delimiter] : delimiters)
      {
        // Loop through all known variableAttributes
        for (auto &[other_index, other_variable] : other_var_attributes)
          {
            // Grab the potential dependency for the current variable. For example, if we
            // are in the variableAttributre for "phi" enumerate all possible dependencies
            // for the variable, like "change(phi)".
            const std::string possible_dependency =
              delimiter.first + other_variable.name + delimiter.second;

            // Populate the dependencies
            if (dependencies.find(possible_dependency) != dependencies.end())
              {
                dependencyType dep_type = relevant_flag.at(variation).first;
                EvalFlags      flags    = relevant_flag.at(variation).second;

                validate_dependency(variation, dep_type, other_index, context);

                std::pair<uint, dependencyType> key = {other_index, dep_type};

                if (eval_flag_set.find(key) != eval_flag_set.end())
                  {
                    eval_flag_set[key] |= flags;
                  }
                else
                  {
                    eval_flag_set.emplace(key, flags);
                  }
              }
          }
      }
  };

  set_dependencies(dependencies_RHS, eval_flag_set_RHS, "RHS");
  set_dependencies(dependencies_LHS, eval_flag_set_LHS, "LHS");

  // Compute the dependency_set
  compute_dependency_set(other_var_attributes);
}

void
variableAttributes::determine_field_solve_type(AttributesList &other_var_attributes,
                                               const bool     &is_postprocessed)
{
  Assert(!eval_flag_set_RHS.empty(),
         dealii::ExcMessage("To begin determining the field solve type, the "
                            "eval_flag_set_RHS must be populated."));

  // Early return for explicit solve types
  if (pde_type == PDEType::EXPLICIT_TIME_DEPENDENT)
    {
      is_postprocessed ? field_solve_type = fieldSolveType::EXPLICIT_POSTPROCESS
                       : field_solve_type = fieldSolveType::EXPLICIT;
      return;
    }

  // Early return for constant fields
  if (pde_type == PDEType::CONSTANT)
    {
      field_solve_type = fieldSolveType::EXPLICIT_CONSTANT;
      return;
    }

  // Determine if we co-nonlinearity. In other words, if the solve of this field is
  // coupled with the solve of one or more field. This takes precedence over the other
  // classifications such as SELF_NONLINEAR.
  //
  // This co-nonlinearity problem is a cycle detection problem. We employ a simple
  // depth-first search algorithm to determine which fields are co-nonlinear. Currently,
  // all co-nonlinear fields are lumped together and solved. When dealing with two sets of
  // co-nonlinear fields performance may suffer greatly.
  find_circular_dependencies(other_var_attributes);
  if (field_solve_type == fieldSolveType::NONEXPLICIT_CO_NONLINEAR)
    {
      return;
    }

  // Check for self-nonlinear solves.
  bool LHS_has_change_eval_flags = false;
  bool LHS_has_normal_eval_flags = false;

  LHS_has_change_eval_flags =
    eval_flag_set_LHS.find({field_index, dependencyType::CHANGE}) !=
    eval_flag_set_LHS.end();
  LHS_has_normal_eval_flags =
    eval_flag_set_LHS.find({field_index, dependencyType::NORMAL}) !=
    eval_flag_set_LHS.end();

  if (LHS_has_change_eval_flags && LHS_has_normal_eval_flags)
    {
      field_solve_type = fieldSolveType::NONEXPLICIT_SELF_NONLINEAR;
      return;
    }

  // Lastly, assign plain linear & auxiliary solves
  if (pde_type == PDEType::AUXILIARY)
    {
      field_solve_type = fieldSolveType::NONEXPLICIT_AUXILIARY;
      return;
    }
  else if (pde_type == PDEType::TIME_INDEPENDENT ||
           pde_type == PDEType::IMPLICIT_TIME_DEPENDENT)
    {
      field_solve_type = fieldSolveType::NONEXPLICIT_LINEAR;
      return;
    }

  // Set undefined is anything falls through our criteria.
  field_solve_type = fieldSolveType::UNDEFINED_SOLVE;
}

void
variableAttributes::print() const
{
  conditionalOStreams::pout_summary()
    << "================================================\n"
    << "\tVariable attribute for " << name << "\n"
    << "================================================\n"
    << "Name: " << name << "\n"
    << "Index: " << field_index << "\n"
    << "Variable type: " << to_string(field_type) << "\n"
    << "Equation type: " << to_string(pde_type) << "\n"
    << "Field solve type: " << to_string(field_solve_type) << "\n";

  conditionalOStreams::pout_summary() << "Evaluation flags RHS:\n";
  for (const auto &[key, value] : eval_flag_set_RHS)
    {
      conditionalOStreams::pout_summary()
        << "  Index: " << key.first << "\n"
        << "  Dependency type: " << to_string(key.second) << "\n"
        << "  Evaluation flags: " << eval_flags_to_string(value) << "\n\n";
    }

  conditionalOStreams::pout_summary() << "Evaluation flags LHS:\n";
  for (const auto &[key, value] : eval_flag_set_LHS)
    {
      conditionalOStreams::pout_summary()
        << "  Index: " << key.first << "\n"
        << "  Dependency type: " << to_string(key.second) << "\n"
        << "  Evaluation flags: " << eval_flags_to_string(value) << "\n\n";
    }

  conditionalOStreams::pout_summary()
    << "Residual flags RHS: " << eval_flags_to_string(eval_flags_residual_RHS) << "\n"
    << "Residual flags LHS: " << eval_flags_to_string(eval_flags_residual_LHS) << "\n"
    << "\n"
    << std::flush;
}

void
variableAttributes::validate_dependency([[maybe_unused]] const std::string &variation,
                                        [[maybe_unused]] dependencyType     dep_type,
                                        [[maybe_unused]] const uint        &other_index,
                                        [[maybe_unused]] const std::string &context) const
{
  Assert(context != "RHS" || dep_type != dependencyType::CHANGE,
         dealii::ExcMessage("Dependencies with the delimiter change(var) are "
                            "not allowed on the RHS of any PDE."));
  Assert(context != "LHS" || pde_type == PDEType::IMPLICIT_TIME_DEPENDENT ||
           pde_type == PDEType::TIME_INDEPENDENT,
         dealii::ExcMessage(
           "Only `TIME_INDEPENDENT` and `IMPLICIT_TIME_DEPENDENT` fields may "
           "have dependencies on the LHS."));
  Assert(dep_type != dependencyType::CHANGE || other_index == field_index,
         dealii::ExcMessage(
           "Dependencies with the delimiter change(var) are only allowed as "
           "dependencies for the same field (e.g, change(phi) is only "
           "allowed as a dependency for phi)."));
  Assert(field_type == fieldType::VECTOR ||
           variation.find("divergence") == std::string::npos,
         dealii::ExcMessage("Dependencies with the divergence delimiter are "
                            "only allowed on vector fields."));
  Assert(field_type == fieldType::VECTOR ||
           variation.find("symmetric_gradient") == std::string::npos,
         dealii::ExcMessage("Dependencies with the symmetric gradient delimiter are "
                            "only allowed on vector fields."));
  Assert(field_type == fieldType::VECTOR || variation.find("curl") == std::string::npos,
         dealii::ExcMessage("Dependencies with the curl delimiter are "
                            "only allowed on vector fields."));
}

void
variableAttributes::compute_dependency_set(const AttributesList &other_var_attributes)
{
  // Compute dependencies for a given eval_flag_set. Change and old flags are irrelevant,
  // so ignore those. If we have EvalFlags::nothing ignore it. If the dependency is
  // explicit or constant ignore it. If the dependency is itself ignore it so it doesn't
  // interfere with the map and flag circularity.
  auto compute_dependencies = [&](const auto &eval_flag_set)
  {
    for (const auto &[pair, flag] : eval_flag_set)
      {
        if (pair.second != dependencyType::NORMAL || flag == EvalFlags::nothing ||
            other_var_attributes.at(pair.first).pde_type ==
              PDEType::EXPLICIT_TIME_DEPENDENT ||
            other_var_attributes.at(pair.first).pde_type == PDEType::CONSTANT ||
            pair.first == field_index)
          {
            continue;
          }

        dependency_set.insert(pair.first);
      }
  };

  compute_dependencies(eval_flag_set_RHS);
  compute_dependencies(eval_flag_set_LHS);
}

void
variableAttributes::find_circular_dependencies(const AttributesList &other_var_attributes)
{
  // Set for visited nodes
  std::set<uint> visited;

  // Set of nodes in the current recursion stack
  std::set<uint> current_stack;

  recursive_DFS(other_var_attributes, visited, current_stack, field_index);
}

void
variableAttributes::recursive_DFS(const AttributesList &other_var_attributes,
                                  std::set<uint>       &visited,
                                  std::set<uint>       &current_stack,
                                  const uint           &vertex)
{
  // Insert the current vertex
  visited.insert(vertex);
  current_stack.insert(vertex);

  // Loop over the dependencies of the current field (the connected nodes)
  for (const auto &dependency : other_var_attributes.at(vertex).dependency_set)
    {
      // If the current recursion stack already has the dependency, we have a cycle
      if (current_stack.find(dependency) != current_stack.end())
        {
          field_solve_type = fieldSolveType::NONEXPLICIT_CO_NONLINEAR;
          return;
        }
      // Otherwise, if we haven't already visited this node continue down the graph
      if (visited.find(dependency) == visited.end())
        {
          recursive_DFS(other_var_attributes, visited, current_stack, dependency);
        }
    }

  // Remove the node from the current recursion stack
  current_stack.erase(vertex);
}
