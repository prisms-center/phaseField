#include "../../include/EquationDependencyParser.h"

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>

#include <iostream>

void
EquationDependencyParser::strip_dependency_whitespace(std::string &dependency_list)
{
  dependency_list.erase(std::remove(dependency_list.begin(), dependency_list.end(), ' '),
                        dependency_list.end());
}

void
EquationDependencyParser::parse(std::vector<std::string> &var_name,
                                std::vector<PDEType>      var_eq_type,
                                std::vector<std::string>  sorted_dependencies_value_RHS,
                                std::vector<std::string> sorted_dependencies_gradient_RHS,
                                std::vector<std::string> sorted_dependencies_value_LHS,
                                std::vector<std::string> sorted_dependencies_gradient_LHS,
                                std::vector<bool>       &var_nonlinear)
{
  // Determine the number of variables
  size_t n_variables = var_name.size();

  // Resize the dependency evaluation flag vectors
  eval_flags_explicit_RHS.resize(n_variables, dealii::EvaluationFlags::nothing);
  eval_flags_nonexplicit_RHS.resize(n_variables, dealii::EvaluationFlags::nothing);
  eval_flags_nonexplicit_LHS.resize(n_variables, dealii::EvaluationFlags::nothing);
  eval_flags_change_nonexplicit_LHS.resize(n_variables, dealii::EvaluationFlags::nothing);

  // Resize the residual evaluation flag vectors
  eval_flags_residual_explicit_RHS.resize(n_variables, dealii::EvaluationFlags::nothing);
  eval_flags_residual_nonexplicit_RHS.resize(n_variables,
                                             dealii::EvaluationFlags::nothing);
  eval_flags_residual_nonexplicit_LHS.resize(n_variables,
                                             dealii::EvaluationFlags::nothing);

  // Now parse the dependency strings to set the flags to true where needed
  for (unsigned int i = 0; i < var_name.size(); i++)
    {
      // Strip excess whitespace
      strip_dependency_whitespace(sorted_dependencies_value_RHS[i]);
      strip_dependency_whitespace(sorted_dependencies_gradient_RHS[i]);

      // Now check for each variable_eq_type
      if (var_eq_type[i] == EXPLICIT_TIME_DEPENDENT)
        {
          bool single_var_nonlinear;

          parseDependencyListRHS(var_name,
                                 var_eq_type,
                                 i,
                                 sorted_dependencies_value_RHS[i],
                                 sorted_dependencies_gradient_RHS[i],
                                 eval_flags_explicit_RHS,
                                 eval_flags_residual_explicit_RHS,
                                 single_var_nonlinear);

          var_nonlinear.push_back(single_var_nonlinear);
        }
      else if (var_eq_type[i] == AUXILIARY)
        {
          bool single_var_nonlinear;

          parseDependencyListRHS(var_name,
                                 var_eq_type,
                                 i,
                                 sorted_dependencies_value_RHS[i],
                                 sorted_dependencies_gradient_RHS[i],
                                 eval_flags_nonexplicit_RHS,
                                 eval_flags_residual_nonexplicit_RHS,
                                 single_var_nonlinear);

          var_nonlinear.push_back(single_var_nonlinear);
        }
      else if (var_eq_type[i] == IMPLICIT_TIME_DEPENDENT ||
               var_eq_type[i] == TIME_INDEPENDENT)
        {
          bool single_var_nonlinear_RHS, single_var_nonlinear_LHS;

          parseDependencyListRHS(var_name,
                                 var_eq_type,
                                 i,
                                 sorted_dependencies_value_RHS[i],
                                 sorted_dependencies_gradient_RHS[i],
                                 eval_flags_nonexplicit_RHS,
                                 eval_flags_residual_nonexplicit_RHS,
                                 single_var_nonlinear_RHS);

          parseDependencyListLHS(var_name,
                                 var_eq_type,
                                 i,
                                 sorted_dependencies_value_LHS[i],
                                 sorted_dependencies_gradient_LHS[i],
                                 eval_flags_nonexplicit_LHS,
                                 eval_flags_change_nonexplicit_LHS,
                                 eval_flags_residual_nonexplicit_LHS,
                                 single_var_nonlinear_LHS);

          var_nonlinear.push_back(single_var_nonlinear_RHS || single_var_nonlinear_LHS);
        }
    }
}

void
EquationDependencyParser::parseDependencyListRHS(
  std::vector<std::string>                              &variable_name_list,
  std::vector<PDEType>                                   variable_eq_type,
  unsigned int                                           variable_index,
  std::string                                           &value_dependencies,
  std::string                                           &gradient_dependencies,
  std::vector<dealii::EvaluationFlags::EvaluationFlags> &evaluation_flags,
  std::vector<dealii::EvaluationFlags::EvaluationFlags> &residual_flags,
  bool                                                  &is_nonlinear)
{
  // Split the dependency strings into lists of entries
  std::vector<std::string> split_value_dependency_list =
    dealii::Utilities::split_string_list(value_dependencies);
  std::vector<std::string> split_gradient_dependency_list =
    dealii::Utilities::split_string_list(gradient_dependencies);

  // Check if either is empty and set value and gradient flags for the
  // residual appropriately
  if (split_value_dependency_list.size() > 0)
    {
      residual_flags[variable_index] |= dealii::EvaluationFlags::values;
    }
  if (split_gradient_dependency_list.size() > 0)
    {
      residual_flags[variable_index] |= dealii::EvaluationFlags::gradients;
    }

  // Merge the lists of dependency entries
  std::vector<std::string> split_dependency_list = split_value_dependency_list;
  split_dependency_list.insert(split_dependency_list.end(),
                               split_gradient_dependency_list.begin(),
                               split_gradient_dependency_list.end());

  // Set nonlinearity to false
  is_nonlinear = false;

  // Cycle through each dependency entry
  for (const auto &dependency : split_dependency_list)
    {
      // Flag to make sure we have assigned a dependency entry
      [[maybe_unused]] bool dependency_entry_assigned = false;

      // Loop through all known variable names [x, grad(x), and hess(x)] to see which ones
      // are on our dependency list. If we have two variables x and y this will loop twice
      // to see if the supplied dependency needs either two of the variables. A successful
      // match will update the values/gradient/hessian flag for that dependency variable.
      std::size_t dependency_variable_index = 0;
      for (const auto &variable : variable_name_list)
        {
          // Create grad() and hess() variants of the variable name
          std::string gradient_variable = {"grad()"};
          gradient_variable.insert(--gradient_variable.end(),
                                   variable.begin(),
                                   variable.end());

          std::string hessian_variable = {"hess()"};
          hessian_variable.insert(--hessian_variable.end(),
                                  variable.begin(),
                                  variable.end());

          // Is the variable we are finding the dependencies for explicit
          bool variable_is_explicit =
            variable_eq_type[variable_index] == EXPLICIT_TIME_DEPENDENT;

          // Is the dependency variable explicit
          bool dependency_variable_is_explicit =
            variable_eq_type[dependency_variable_index] == EXPLICIT_TIME_DEPENDENT;

          // Is the dependency the variable
          bool same_variable = variable_index == dependency_variable_index;

          // Case if the dependency is x
          if (dependency == variable)
            {
              evaluation_flags[dependency_variable_index] |=
                dealii::EvaluationFlags::values;
              dependency_entry_assigned = true;
            }
          // Case if the dependency is grad(x)
          else if (dependency == gradient_variable)
            {
              evaluation_flags[dependency_variable_index] |=
                dealii::EvaluationFlags::gradients;
              dependency_entry_assigned = true;
            }
          // Case if the dependency is hess(x)
          else if (dependency == hessian_variable)
            {
              evaluation_flags[dependency_variable_index] |=
                dealii::EvaluationFlags::hessians;
              dependency_entry_assigned = true;
            }

          // Check for nonlinearity
          is_nonlinear =
            !variable_is_explicit && !same_variable && !dependency_variable_is_explicit;

          // Increment counter
          ++dependency_variable_index;
        }

      Assert(dependency_entry_assigned,
             dealii::StandardExceptions::ExcMessage("PRISMS-PF Error: Dependency entry " +
                                                    dependency + " is not valid."));
    }
}

void
EquationDependencyParser::parseDependencyListLHS(
  std::vector<std::string>                              &variable_name_list,
  std::vector<PDEType>                                   variable_eq_type,
  unsigned int                                           variable_index,
  std::string                                           &value_dependencies,
  std::string                                           &gradient_dependencies,
  std::vector<dealii::EvaluationFlags::EvaluationFlags> &evaluation_flags,
  std::vector<dealii::EvaluationFlags::EvaluationFlags> &change_flags,
  std::vector<dealii::EvaluationFlags::EvaluationFlags> &residual_flags,
  bool                                                  &is_nonlinear)
{
  // Split the dependency strings into lists of entries
  std::vector<std::string> split_value_dependency_list =
    dealii::Utilities::split_string_list(value_dependencies);
  std::vector<std::string> split_gradient_dependency_list =
    dealii::Utilities::split_string_list(gradient_dependencies);

  // Check if either is empty and set value and gradient flags for the
  // residual appropriately
  if (split_value_dependency_list.size() > 0)
    {
      residual_flags[variable_index] |= dealii::EvaluationFlags::values;
    }
  if (split_gradient_dependency_list.size() > 0)
    {
      residual_flags[variable_index] |= dealii::EvaluationFlags::gradients;
    }

  // Merge the lists of dependency entries
  std::vector<std::string> split_dependency_list = split_value_dependency_list;
  split_dependency_list.insert(split_dependency_list.end(),
                               split_gradient_dependency_list.begin(),
                               split_gradient_dependency_list.end());

  // Set nonlinearity to false
  is_nonlinear = false;

  // Cycle through each dependency entry
  for (const auto &dependency : split_dependency_list)
    {
      // Flag to make sure we have assigned a dependency entry
      [[maybe_unused]] bool dependency_entry_assigned = false;

      // Loop through all known variable names [x, grad(x), and hess(x)] to see which ones
      // are on our dependency list. If we have two variables x and y this will loop twice
      // to see if the supplied dependency needs either two of the variables. A successful
      // match will update the values/gradient/hessian flag for that dependency variable.
      std::size_t dependency_variable_index = 0;
      for (const auto &variable : variable_name_list)
        {
          // Create grad(), hess(), change(), grad(change()), and hess(change()) variants
          // of the variable name
          std::string gradient_variable = {"grad()"};
          gradient_variable.insert(--gradient_variable.end(),
                                   variable.begin(),
                                   variable.end());

          std::string hessian_variable = {"hess()"};
          hessian_variable.insert(--hessian_variable.end(),
                                  variable.begin(),
                                  variable.end());

          std::string change_value_variable = {"change()"};
          change_value_variable.insert(--change_value_variable.end(),
                                       variable.begin(),
                                       variable.end());

          std::string change_gradient_variable = {"grad(change())"};
          change_gradient_variable.insert(--(--change_gradient_variable.end()),
                                          variable.begin(),
                                          variable.end());

          std::string change_hessian_variable = {"hess(change())"};
          change_hessian_variable.insert(--(--change_hessian_variable.end()),
                                         variable.begin(),
                                         variable.end());

          // Is the variable we are finding the dependencies for explicit
          bool dependency_variable_is_explicit =
            variable_eq_type[dependency_variable_index] == EXPLICIT_TIME_DEPENDENT;

          // Case if the dependency is x
          if (dependency == variable)
            {
              evaluation_flags[dependency_variable_index] |=
                dealii::EvaluationFlags::values;
              dependency_entry_assigned = true;

              // Check for nonlinearity
              is_nonlinear = !dependency_variable_is_explicit;
            }
          // Case if the dependency is grad(x)
          else if (dependency == gradient_variable)
            {
              evaluation_flags[dependency_variable_index] |=
                dealii::EvaluationFlags::gradients;
              dependency_entry_assigned = true;

              // Check for nonlinearity
              is_nonlinear = !dependency_variable_is_explicit;
            }
          // Case if the dependency is hess(x)
          else if (dependency == hessian_variable)
            {
              evaluation_flags[dependency_variable_index] |=
                dealii::EvaluationFlags::hessians;
              dependency_entry_assigned = true;

              // Check for nonlinearity
              is_nonlinear = !dependency_variable_is_explicit;
            }
          // Case if the dependency is change(x)
          else if (dependency == change_value_variable)
            {
              change_flags[dependency_variable_index] |= dealii::EvaluationFlags::values;
              dependency_entry_assigned = true;

              Assert(variable_index == dependency_variable_index,
                     dealii::StandardExceptions::ExcMessage(
                       "PRISMS-PF Error: Dependency entry " + dependency +
                       " is not valid because the change in a variable can "
                       "only be accessed in its own governing equation."));
            }
          // Case if the dependency is grad(change(x))
          else if (dependency == change_gradient_variable)
            {
              change_flags[dependency_variable_index] |=
                dealii::EvaluationFlags::gradients;
              dependency_entry_assigned = true;

              Assert(variable_index == dependency_variable_index,
                     dealii::StandardExceptions::ExcMessage(
                       "PRISMS-PF Error: Dependency entry " + dependency +
                       " is not valid because the change in a variable can "
                       "only be accessed in its own governing equation."));
            }
          // Case if the dependency is hess(change(x))
          else if (dependency == change_hessian_variable)
            {
              change_flags[dependency_variable_index] |=
                dealii::EvaluationFlags::hessians;
              dependency_entry_assigned = true;

              Assert(variable_index == dependency_variable_index,
                     dealii::StandardExceptions::ExcMessage(
                       "PRISMS-PF Error: Dependency entry " + dependency +
                       " is not valid because the change in a variable can "
                       "only be accessed in its own governing equation."));
            }

          // Increment counter
          ++dependency_variable_index;
        }

      Assert(dependency_entry_assigned,
             dealii::StandardExceptions::ExcMessage("PRISMS-PF Error: Dependency entry " +
                                                    dependency + " is not valid."));
    }
}

void
EquationDependencyParser::pp_parse(std::vector<std::string> &var_name,
                                   std::vector<std::string> &pp_var_name,
                                   std::vector<std::string>  sorted_dependencies_value,
                                   std::vector<std::string>  sorted_dependencies_gradient)
{
  // Determine the number of variables
  size_t n_variables             = var_name.size();
  size_t n_postprocess_variables = pp_var_name.size();

  // Resize the dependency evaluation flag vectors
  eval_flags_postprocess.resize(n_variables, dealii::EvaluationFlags::nothing);

  // Resize the residual evaluation flag vectors
  eval_flags_residual_postprocess.resize(n_postprocess_variables,
                                         dealii::EvaluationFlags::nothing);

  // Now parse the dependency strings to set the flags to true where needed
  for (unsigned int i = 0; i < pp_var_name.size(); i++)
    {
      // Strip excess whitespace
      strip_dependency_whitespace(sorted_dependencies_value[i]);
      strip_dependency_whitespace(sorted_dependencies_gradient[i]);

      parseDependencyListPP(var_name,
                            i,
                            sorted_dependencies_value[i],
                            sorted_dependencies_gradient[i],
                            eval_flags_postprocess,
                            eval_flags_residual_postprocess);
    }
}

void
EquationDependencyParser::parseDependencyListPP(
  std::vector<std::string>                              &variable_name_list,
  unsigned int                                           variable_index,
  std::string                                           &value_dependencies,
  std::string                                           &gradient_dependencies,
  std::vector<dealii::EvaluationFlags::EvaluationFlags> &evaluation_flags,
  std::vector<dealii::EvaluationFlags::EvaluationFlags> &residual_flags)
{
  // Split the dependency strings into lists of entries
  std::vector<std::string> split_value_dependency_list =
    dealii::Utilities::split_string_list(value_dependencies);
  std::vector<std::string> split_gradient_dependency_list =
    dealii::Utilities::split_string_list(gradient_dependencies);

  // Check if either is empty and set value and gradient flags for the
  // residual appropriately
  if (split_value_dependency_list.size() > 0)
    {
      residual_flags[variable_index] |= dealii::EvaluationFlags::values;
    }
  if (split_gradient_dependency_list.size() > 0)
    {
      residual_flags[variable_index] |= dealii::EvaluationFlags::gradients;
    }

  // Merge the lists of dependency entries
  std::vector<std::string> split_dependency_list = split_value_dependency_list;
  split_dependency_list.insert(split_dependency_list.end(),
                               split_gradient_dependency_list.begin(),
                               split_gradient_dependency_list.end());

  // Cycle through each dependency entry
  for (const auto &dependency : split_dependency_list)
    {
      // Flag to make sure we have assigned a dependency entry
      [[maybe_unused]] bool dependency_entry_assigned = false;

      // Loop through all known variable names [x, grad(x), and hess(x)] to see which ones
      // are on our dependency list. If we have two variables x and y this will loop twice
      // to see if the supplied dependency needs either two of the variables. A successful
      // match will update the values/gradient/hessian flag for that dependency variable.
      std::size_t dependency_variable_index = 0;
      for (const auto &variable : variable_name_list)
        {
          // Create grad() and hess() variants of the variable name
          std::string gradient_variable = {"grad()"};
          gradient_variable.insert(--gradient_variable.end(),
                                   variable.begin(),
                                   variable.end());

          std::string hessian_variable = {"hess()"};
          hessian_variable.insert(--hessian_variable.end(),
                                  variable.begin(),
                                  variable.end());

          if (dependency == variable)
            {
              evaluation_flags[dependency_variable_index] |=
                dealii::EvaluationFlags::values;
              dependency_entry_assigned = true;
            }
          else if (dependency == gradient_variable)
            {
              evaluation_flags[dependency_variable_index] |=
                dealii::EvaluationFlags::gradients;
              dependency_entry_assigned = true;
            }
          else if (dependency == hessian_variable)
            {
              evaluation_flags[dependency_variable_index] |=
                dealii::EvaluationFlags::hessians;
              dependency_entry_assigned = true;
            }

          // Increment counter
          ++dependency_variable_index;
        }

      Assert(dependency_entry_assigned,
             dealii::StandardExceptions::ExcMessage("PRISMS-PF Error: Dependency entry " +
                                                    dependency + " is not valid."));
    }
}
