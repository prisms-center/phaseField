// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include <map>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

void
VariableAttributes::format_dependencies()
{
  raw_dependencies.dependencies_rhs.insert(
    raw_dependencies.dependencies_value_rhs.begin(),
    raw_dependencies.dependencies_value_rhs.end());
  raw_dependencies.dependencies_rhs.insert(
    raw_dependencies.dependencies_gradient_rhs.begin(),
    raw_dependencies.dependencies_gradient_rhs.end());

  raw_dependencies.dependencies_lhs.insert(
    raw_dependencies.dependencies_value_lhs.begin(),
    raw_dependencies.dependencies_value_lhs.end());
  raw_dependencies.dependencies_lhs.insert(
    raw_dependencies.dependencies_gradient_lhs.begin(),
    raw_dependencies.dependencies_gradient_lhs.end());
}

void
VariableAttributes::parse_residual_dependencies()
{
  // Check if either is empty and set value and gradient flags for the
  // residual appropriately. Note that this relies on the fact that only valid
  // dependencies are part of the set.
  if (!raw_dependencies.dependencies_value_rhs.empty())
    {
      eval_flags_residual_rhs |= dealii::EvaluationFlags::EvaluationFlags::values;
    }
  if (!raw_dependencies.dependencies_gradient_rhs.empty())
    {
      eval_flags_residual_rhs |= dealii::EvaluationFlags::EvaluationFlags::gradients;
    }
  if (!raw_dependencies.dependencies_value_lhs.empty())
    {
      eval_flags_residual_lhs |= dealii::EvaluationFlags::EvaluationFlags::values;
    }
  if (!raw_dependencies.dependencies_gradient_lhs.empty())
    {
      eval_flags_residual_lhs |= dealii::EvaluationFlags::EvaluationFlags::gradients;
    }
}

void
VariableAttributes::parse_dependencies(
  std::map<unsigned int, VariableAttributes> &other_var_attributes,
  const Types::Index                         &_max_fields,
  const Types::Index                         &_max_dependency_types)
{
  const std::map<std::string, std::pair<DependencyType, EvalFlags>> relevant_flag = []()
  {
    // Modifiers for the dependency types
    const std::vector<std::pair<std::string, DependencyType>> dependency_types = {
      {"",        DependencyType::Normal  },
      {"_change", DependencyType::Change  },
      {"_old_1",  DependencyType::OldOne  },
      {"_old_2",  DependencyType::OldTwo  },
      {"_old_3",  DependencyType::OldThree},
      {"_old_4",  DependencyType::OldFour },
    };

    // Dependency evaluation types & their corresponding evaluation flag
    const std::vector<std::pair<std::string, EvalFlags>> eval_flags = {
      {"value",              dealii::EvaluationFlags::EvaluationFlags::values   },
      {"gradient",           dealii::EvaluationFlags::EvaluationFlags::gradients},
      {"hessian",            dealii::EvaluationFlags::EvaluationFlags::hessians },
      {"hessian_diagonal",   dealii::EvaluationFlags::EvaluationFlags::hessians },
      {"laplacian",          dealii::EvaluationFlags::EvaluationFlags::hessians },
      {"divergence",         dealii::EvaluationFlags::EvaluationFlags::gradients},
      {"symmetric_gradient", dealii::EvaluationFlags::EvaluationFlags::gradients},
      {"curl",               dealii::EvaluationFlags::EvaluationFlags::gradients},
    };

    std::map<std::string, std::pair<DependencyType, EvalFlags>> map;
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
    const std::vector<std::pair<std::string, std::pair<std::string, std::string>>>
      dependency_types_delimiters = {
        {"",        {"", ""}        },
        {"_change", {"change(", ")"}},
        {"_old_1",  {"old_1(", ")"} },
        {"_old_2",  {"old_2(", ")"} },
        {"_old_3",  {"old_3(", ")"} },
        {"_old_4",  {"old_4(", ")"} },
    };

    // Dependency evaluation types & their corresponding delimiter.
    const std::vector<std::pair<std::string, std::pair<std::string, std::string>>>
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
  auto set_dependencies = [&](const std::set<std::string>         &dependencies,
                              std::vector<std::vector<EvalFlags>> &eval_flag_set,
                              const std::string                   &context)
  {
    // Loop through the available delimiters
    for (const auto &[variation, delimiter] : delimiters)
      {
        // Loop through all known VariableAttributes
        for (auto &[other_index, other_variable] : other_var_attributes)
          {
            // Grab the potential dependency for the variable. For example, if we are in
            // the variableAttributre for "phi" enumerate all possible dependencies for
            // the variable, like "change(phi)".
            const std::string possible_dependency =
              delimiter.first + other_variable.name + delimiter.second;

            // Populate the dependencies
            if (dependencies.contains(possible_dependency))
              {
                const DependencyType        dep_type = relevant_flag.at(variation).first;
                const EvalFlags             flags    = relevant_flag.at(variation).second;
                const FieldInfo::TensorRank other_field_type =
                  other_variable.field_info.tensor_rank;

                validate_dependency(variation,
                                    dep_type,
                                    other_index,
                                    other_field_type,
                                    context);

                eval_flag_set[other_index][static_cast<Types::Index>(dep_type)] |= flags;
              }
          }
      }
  };

  // Initialize the eval flag sets
  eval_flag_set_rhs.resize(
    _max_fields,
    std::vector<EvalFlags>(_max_dependency_types + 1,
                           dealii::EvaluationFlags::EvaluationFlags::nothing));
  eval_flag_set_lhs.resize(
    _max_fields,
    std::vector<EvalFlags>(_max_dependency_types + 1,
                           dealii::EvaluationFlags::EvaluationFlags::nothing));

  set_dependencies(raw_dependencies.dependencies_rhs, eval_flag_set_rhs, "RHS");
  set_dependencies(raw_dependencies.dependencies_lhs, eval_flag_set_lhs, "LHS");

  for (const auto &[other_index, other_variable] : other_var_attributes)
    {
      if (raw_dependencies.nucleating_fields.contains(other_variable.name))
        {
          nucleation_indices.push_back(other_index);
        }
    }

  // Compute the dependency_set and simplified_dependency_set
  compute_dependency_set(other_var_attributes);
  compute_simplified_dependency_set(other_var_attributes);
}

void
VariableAttributes::determine_field_solve_type(
  const std::map<unsigned int, VariableAttributes> &other_var_attributes)
{
  // Early return for constant fields
  if (pde_type == PDEType::Constant)
    {
      field_solve_type = FieldSolveType::ExplicitConstant;
      return;
    }

  AssertThrow(!eval_flag_set_rhs.empty(),
              dealii::ExcMessage(
                "To begin determining the field solve type that is not constant, the "
                "eval_flag_set_rhs must be populated."));

  // Early return for explicit solve types
  if (pde_type == PDEType::ExplicitTimeDependent)
    {
      is_postprocessed_variable ? field_solve_type = FieldSolveType::ExplicitPostprocess
                                : field_solve_type = FieldSolveType::Explicit;
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
  if (field_solve_type == FieldSolveType::NonexplicitCononlinear)
    {
      return;
    }

  // Check for self-nonlinear solves.
  if (eval_flag_set_lhs[field_index][static_cast<Types::Index>(DependencyType::Change)] !=
        EvalFlags::nothing &&
      eval_flag_set_lhs[field_index][static_cast<Types::Index>(DependencyType::Normal)] !=
        EvalFlags::nothing)
    {
      field_solve_type = FieldSolveType::NonexplicitSelfnonlinear;
      return;
    }

  // Lastly, assign plain linear & auxiliary solves
  if (pde_type == PDEType::Auxiliary)
    {
      field_solve_type = FieldSolveType::NonexplicitAuxiliary;
      return;
    }
  if (pde_type == PDEType::TimeIndependent || pde_type == PDEType::ImplicitTimeDependent)
    {
      field_solve_type = FieldSolveType::NonexplicitLinear;
      return;
    }

  // Set undefined is anything falls through our criteria.
  field_solve_type = Numbers::invalid_field_solve_type;
}

void
VariableAttributes::print() const
{
  ConditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Variable attribute for " << name << "\n"
    << "================================================\n"
    << "Name: " << name << "\n"
    << "Index: " << field_index << "\n"
    << "Variable type: " << to_string(field_info.tensor_rank) << "\n"
    << "Equation type: " << to_string(pde_type) << "\n"
    << "Postprocessed field: " << bool_to_string(is_postprocessed_variable) << "\n"
    << "Field solve type: " << to_string(field_solve_type) << "\n";

  ConditionalOStreams::pout_summary() << "Evaluation flags RHS:\n";
  Types::Index index = 0;
  for (const auto &dependency_set : eval_flag_set_rhs)
    {
      Types::Index dep_index = 0;
      for (const auto &value : dependency_set)
        {
          // Skip where the evaluation flags are nothing
          if (value == dealii::EvaluationFlags::EvaluationFlags::nothing)
            {
              dep_index++;
              continue;
            }

          ConditionalOStreams::pout_summary()
            << "  Index: " << index << "\n"
            << "  Dependency type: " << to_string(static_cast<DependencyType>(dep_index))
            << "\n"
            << "  Evaluation flags: " << eval_flags_to_string(value) << "\n\n";
          dep_index++;
        }
      index++;
    }

  ConditionalOStreams::pout_summary() << "Evaluation flags LHS:\n";
  index = 0;
  for (const auto &dependency_set : eval_flag_set_lhs)
    {
      Types::Index dep_index = 0;
      for (const auto &value : dependency_set)
        {
          // Skip where the evaluation flags are nothing
          if (value == dealii::EvaluationFlags::EvaluationFlags::nothing)
            {
              dep_index++;
              continue;
            }

          ConditionalOStreams::pout_summary()
            << "  Index: " << index << "\n"
            << "  Dependency type: " << to_string(static_cast<DependencyType>(dep_index))
            << "\n"
            << "  Evaluation flags: " << eval_flags_to_string(value) << "\n\n";
          dep_index++;
        }
      index++;
    }

  ConditionalOStreams::pout_summary()
    << "Residual flags RHS: " << eval_flags_to_string(eval_flags_residual_rhs) << "\n"
    << "Residual flags LHS: " << eval_flags_to_string(eval_flags_residual_lhs) << "\n"
    << "\n"
    << std::flush;
}

void
VariableAttributes::validate_dependency(
  [[maybe_unused]] const std::string           &variation,
  [[maybe_unused]] DependencyType               dep_type,
  [[maybe_unused]] const unsigned int          &other_index,
  [[maybe_unused]] const FieldInfo::TensorRank &other_field_type,
  [[maybe_unused]] const std::string           &context) const
{
  AssertThrow(context != "RHS" || dep_type != DependencyType::Change,
              dealii::ExcMessage("Dependencies with the delimiter change(var) are "
                                 "not allowed on the RHS of any PDE."));
  AssertThrow(context != "LHS" || pde_type == PDEType::ImplicitTimeDependent ||
                pde_type == PDEType::TimeIndependent,
              dealii::ExcMessage(
                "Only `TimeIndependent` and `ImplicitTimeDependent` fields may "
                "have dependencies on the LHS."));
  AssertThrow(dep_type != DependencyType::Change || other_index == field_index,
              dealii::ExcMessage(
                "Dependencies with the delimiter change(var) are only allowed as "
                "dependencies for the same field (e.g, change(phi) is only "
                "allowed as a dependency for phi)."));
  AssertThrow(other_field_type == FieldInfo::TensorRank::Vector ||
                variation.find("divergence") == std::string::npos,
              dealii::ExcMessage("Dependencies with the divergence delimiter are "
                                 "only allowed on vector fields."));
  AssertThrow(other_field_type == FieldInfo::TensorRank::Vector ||
                variation.find("symmetric_gradient") == std::string::npos,
              dealii::ExcMessage("Dependencies with the symmetric gradient delimiter are "
                                 "only allowed on vector fields."));
  AssertThrow(other_field_type == FieldInfo::TensorRank::Vector ||
                variation.find("curl") == std::string::npos,
              dealii::ExcMessage("Dependencies with the curl delimiter are "
                                 "only allowed on vector fields."));
}

void
VariableAttributes::compute_dependency_set(
  const std::map<unsigned int, VariableAttributes> &other_var_attributes)
{
  // Compute dependencies for a given eval_flag_set. Change flags are irrelevant for RHS,
  // so ignore those. If we have EvalFlags::nothing ignore it. Always add itself as a
  // dependency for RHS since this is used for determining what FEEvaluation objects
  // should be made.

  // First resize the dependency_set_rhs and dependency_set_lhs and populate them with
  // invalid entries. This is so we don't create the FEEvaluation objects.
  dependency_set_rhs.resize(
    eval_flag_set_rhs.size(),
    std::vector<FieldInfo::TensorRank>(eval_flag_set_rhs.begin()->size(),
                                       FieldInfo::TensorRank::Undefined));
  dependency_set_lhs.resize(
    eval_flag_set_lhs.size(),
    std::vector<FieldInfo::TensorRank>(eval_flag_set_lhs.begin()->size(),
                                       FieldInfo::TensorRank::Undefined));

  Types::Index index = 0;
  for (const auto &dependency_set : eval_flag_set_rhs)
    {
      Types::Index dep_index = 0;
      for (const auto &value : dependency_set)
        {
          // TODO (landinjm): Should the change terms be disallowed? The assertion might
          // be redundant but provide some context.
          if (static_cast<DependencyType>(dep_index) == DependencyType::Change ||
              value == EvalFlags::nothing)
            {
              dep_index++;
              continue;
            }

          Assert(other_var_attributes.contains(index),
                 dealii::ExcMessage(
                   "The provided attributes does not have an entry for the index = " +
                   std::to_string(index)));

          dependency_set_rhs[index][dep_index] =
            other_var_attributes.at(index).field_info.tensor_rank;

          dep_index++;
        }

      index++;
    }

  dependency_set_rhs[field_index][static_cast<Types::Index>(DependencyType::Normal)] =
    field_info.tensor_rank;

  index = 0;
  for (const auto &dependency_set : eval_flag_set_lhs)
    {
      Types::Index dep_index = 0;
      for (const auto &value : dependency_set)
        {
          if (value == EvalFlags::nothing)
            {
              dep_index++;
              continue;
            }

          Assert(other_var_attributes.contains(index),
                 dealii::ExcMessage(
                   "The provided attributes does not have an entry for the index = " +
                   std::to_string(index)));

          dependency_set_lhs[index][dep_index] =
            other_var_attributes.at(index).field_info.tensor_rank;

          dep_index++;
        }

      index++;
    }
}

void
VariableAttributes::compute_simplified_dependency_set(
  const std::map<unsigned int, VariableAttributes> &other_var_attributes)
{
  // Compute dependencies for a given eval_flag_set. Change and old flags are irrelevant,
  // so ignore those. If we have EvalFlags::nothing ignore it. If the dependency is
  // explicit or constant ignore it. If the dependency is itself ignore it so it doesn't
  // interfere with the map and flag circularity. If the dependencys have different solve
  // blocks ignore those as well.
  auto compute_dependencies = [&](const auto &eval_flag_set)
  {
    Types::Index index = 0;
    for (const auto &dependency_set : eval_flag_set)
      {
        Types::Index dep_index = 0;
        for (const auto &value : dependency_set)
          {
            if (static_cast<DependencyType>(dep_index) != DependencyType::Normal ||
                value == EvalFlags::nothing ||
                other_var_attributes.at(index).pde_type ==
                  PDEType::ExplicitTimeDependent ||
                other_var_attributes.at(index).pde_type == PDEType::Constant ||
                index == field_index ||
                other_var_attributes.at(index).solve_block !=
                  other_var_attributes.at(dep_index).solve_block)
              {
                dep_index++;
                continue;
              }

            simplified_dependency_set.insert(index);

            dep_index++;
          }

        index++;
      }
  };

  compute_dependencies(eval_flag_set_rhs);
  compute_dependencies(eval_flag_set_lhs);
}

void
VariableAttributes::find_circular_dependencies(
  const std::map<unsigned int, VariableAttributes> &other_var_attributes)
{
  // Set for visited nodes
  std::set<unsigned int> visited;

  // Set of nodes in the current recursion stack
  std::set<unsigned int> current_stack;

  recursive_depth_first_search(other_var_attributes, visited, current_stack, field_index);
}

// NOLINTBEGIN(misc-no-recursion)

void
VariableAttributes::recursive_depth_first_search(
  const std::map<unsigned int, VariableAttributes> &other_var_attributes,
  std::set<unsigned int>                           &visited,
  std::set<unsigned int>                           &current_stack,
  const unsigned int                               &vertex)
{
  // Insert the current vertex
  visited.insert(vertex);
  current_stack.insert(vertex);

  // Loop over the dependencies of the current field (the connected nodes)
  for (const auto &dependency : other_var_attributes.at(vertex).simplified_dependency_set)
    {
      // If the current recursion stack already has the dependency, we have a cycle
      if (current_stack.contains(dependency))
        {
          field_solve_type = FieldSolveType::NonexplicitCononlinear;
          return;
        }
      // Otherwise, if we haven't already visited this node continue down the graph
      if (!visited.contains(dependency))
        {
          recursive_depth_first_search(other_var_attributes,
                                       visited,
                                       current_stack,
                                       dependency);
        }
    }

  // Remove the node from the current recursion stack
  current_stack.erase(vertex);
}

// NOLINTEND(misc-no-recursion)

PRISMS_PF_END_NAMESPACE
