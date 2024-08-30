#include "../../include/EquationDependencyParser.h"

#include <deal.II/base/utilities.h>

#include <iostream>

void
EquationDependencyParser::strip_dependency_whitespace(std::string &dependency_list)
{
  dependency_list.erase(std::remove(dependency_list.begin(), dependency_list.end(), ' '),
                        dependency_list.end());
}

void
EquationDependencyParser::parse(std::vector<std::string> var_name,
                                std::vector<PDEType>     var_eq_type,
                                std::vector<std::string> sorted_dependencies_value_RHS,
                                std::vector<std::string> sorted_dependencies_gradient_RHS,
                                std::vector<std::string> sorted_dependencies_value_LHS,
                                std::vector<std::string> sorted_dependencies_gradient_LHS,
                                std::vector<bool>       &var_nonlinear)
{
  // Determine the number of variables
  size_t n_variables = var_name.size();

  // Resize the dependency vectors
  need_value_explicit_RHS.resize(n_variables, false);
  need_gradient_explicit_RHS.resize(n_variables, false);
  need_hessian_explicit_RHS.resize(n_variables, false);

  need_value_nonexplicit_RHS.resize(n_variables, false);
  need_gradient_nonexplicit_RHS.resize(n_variables, false);
  need_hessian_nonexplicit_RHS.resize(n_variables, false);

  need_value_nonexplicit_LHS.resize(n_variables, false);
  need_gradient_nonexplicit_LHS.resize(n_variables, false);
  need_hessian_nonexplicit_LHS.resize(n_variables, false);

  need_value_change_nonexplicit_LHS.resize(n_variables, false);
  need_gradient_change_nonexplicit_LHS.resize(n_variables, false);
  need_hessian_change_nonexplicit_LHS.resize(n_variables, false);

  // Resize the residual vectors
  need_value_residual_explicit_RHS.resize(n_variables, false);
  need_gradient_residual_explicit_RHS.resize(n_variables, false);

  need_value_residual_nonexplicit_RHS.resize(n_variables, false);
  need_gradient_residual_nonexplicit_RHS.resize(n_variables, false);

  need_value_residual_nonexplicit_LHS.resize(n_variables, false);
  need_gradient_residual_nonexplicit_LHS.resize(n_variables, false);

  // Now parse the dependency strings to set the flags to true where needed
  for (unsigned int i = 0; i < var_name.size(); i++)
    {
      // Strip excess whitespace
      strip_dependency_whitespace(sorted_dependencies_value_RHS[i]);
      strip_dependency_whitespace(sorted_dependencies_gradient_RHS[i]);

      // Now check for each variable_eq_type
      if (var_eq_type[i] == EXPLICIT_TIME_DEPENDENT)
        {
          bool need_value_residual_entry, need_gradient_residual_entry,
            single_var_nonlinear;

          parseDependencyListRHS(var_name,
                                 var_eq_type,
                                 i,
                                 sorted_dependencies_value_RHS[i],
                                 sorted_dependencies_gradient_RHS[i],
                                 need_value_explicit_RHS,
                                 need_gradient_explicit_RHS,
                                 need_hessian_explicit_RHS,
                                 need_value_residual_entry,
                                 need_gradient_residual_entry,
                                 single_var_nonlinear);

          var_nonlinear.push_back(single_var_nonlinear);

          need_value_residual_explicit_RHS[i]    = need_value_residual_entry;
          need_gradient_residual_explicit_RHS[i] = need_gradient_residual_entry;
        }
      else if (var_eq_type[i] == AUXILIARY)
        {
          bool need_value_residual_entry, need_gradient_residual_entry,
            single_var_nonlinear;

          parseDependencyListRHS(var_name,
                                 var_eq_type,
                                 i,
                                 sorted_dependencies_value_RHS[i],
                                 sorted_dependencies_gradient_RHS[i],
                                 need_value_nonexplicit_RHS,
                                 need_gradient_nonexplicit_RHS,
                                 need_hessian_nonexplicit_RHS,
                                 need_value_residual_entry,
                                 need_gradient_residual_entry,
                                 single_var_nonlinear);

          var_nonlinear.push_back(single_var_nonlinear);

          need_value_residual_nonexplicit_RHS[i]    = need_value_residual_entry;
          need_gradient_residual_nonexplicit_RHS[i] = need_gradient_residual_entry;
        }
      else if (var_eq_type[i] == IMPLICIT_TIME_DEPENDENT ||
               var_eq_type[i] == TIME_INDEPENDENT)
        {
          bool need_value_residual_entry, need_gradient_residual_entry,
            single_var_nonlinear_RHS, single_var_nonlinear_LHS;

          parseDependencyListRHS(var_name,
                                 var_eq_type,
                                 i,
                                 sorted_dependencies_value_RHS[i],
                                 sorted_dependencies_gradient_RHS[i],
                                 need_value_nonexplicit_RHS,
                                 need_gradient_nonexplicit_RHS,
                                 need_hessian_nonexplicit_RHS,
                                 need_value_residual_entry,
                                 need_gradient_residual_entry,
                                 single_var_nonlinear_RHS);

          need_value_residual_nonexplicit_RHS[i]    = need_value_residual_entry;
          need_gradient_residual_nonexplicit_RHS[i] = need_gradient_residual_entry;

          parseDependencyListLHS(var_name,
                                 var_eq_type,
                                 i,
                                 sorted_dependencies_value_LHS[i],
                                 sorted_dependencies_gradient_LHS[i],
                                 need_value_nonexplicit_LHS,
                                 need_gradient_nonexplicit_LHS,
                                 need_hessian_nonexplicit_LHS,
                                 need_value_change_nonexplicit_LHS,
                                 need_gradient_change_nonexplicit_LHS,
                                 need_hessian_change_nonexplicit_LHS,
                                 need_value_residual_entry,
                                 need_gradient_residual_entry,
                                 single_var_nonlinear_LHS);

          var_nonlinear.push_back(single_var_nonlinear_RHS || single_var_nonlinear_LHS);

          need_value_residual_nonexplicit_LHS[i]    = need_value_residual_entry;
          need_gradient_residual_nonexplicit_LHS[i] = need_gradient_residual_entry;
        }
    }
}

void
EquationDependencyParser::parseDependencyListRHS(
  std::vector<std::string> variable_name_list,
  std::vector<PDEType>     variable_eq_type,
  unsigned int             variable_index,
  std::string              value_dependencies,
  std::string              gradient_dependencies,
  std::vector<bool>       &need_value,
  std::vector<bool>       &need_gradient,
  std::vector<bool>       &need_hessian,
  bool                    &need_value_residual,
  bool                    &need_gradient_residual,
  bool                    &is_nonlinear)
{
  // Split the dependency strings into lists of entries
  std::vector<std::string> split_value_dependency_list =
    dealii::Utilities::split_string_list(value_dependencies);
  std::vector<std::string> split_gradient_dependency_list =
    dealii::Utilities::split_string_list(gradient_dependencies);

  // Check if either is empty and set need_value_residual and need_gradient
  // residual appropriately
  split_value_dependency_list.size() > 0 ? need_value_residual = true
                                         : need_value_residual = false;
  split_gradient_dependency_list.size() > 0 ? need_gradient_residual = true
                                            : need_gradient_residual = false;

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
      bool dependency_entry_assigned = false;

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
              need_value[dependency_variable_index] = true;
              dependency_entry_assigned             = true;
            }
          // Case if the dependency is grad(x)
          else if (dependency == gradient_variable)
            {
              need_gradient[dependency_variable_index] = true;
              dependency_entry_assigned                = true;
            }
          // Case if the dependency is hess(x)
          else if (dependency == hessian_variable)
            {
              need_hessian.at(dependency_variable_index) = true;
              dependency_entry_assigned                  = true;
            }

          // Check for nonlinearity
          is_nonlinear = is_nonlinear || !variable_is_explicit && !same_variable &&
                                           !dependency_variable_is_explicit;

          // Increment counter
          ++dependency_variable_index;
        }

      Assert(dependency_entry_assigned,
             dealii::StandardExceptions::ExcMessage("PRISMS-PF Error: Dependency entry " +
                                                    dependency + " is not valid."))
    }
}

void
EquationDependencyParser::parseDependencyListLHS(
  std::vector<std::string> variable_name_list,
  std::vector<PDEType>     variable_eq_type,
  unsigned int             variable_index,
  std::string              value_dependencies,
  std::string              gradient_dependencies,
  std::vector<bool>       &need_value,
  std::vector<bool>       &need_gradient,
  std::vector<bool>       &need_hessian,
  std::vector<bool>       &need_value_change,
  std::vector<bool>       &need_gradient_change,
  std::vector<bool>       &need_hessian_change,
  bool                    &need_value_residual,
  bool                    &need_gradient_residual,
  bool                    &is_nonlinear)
{
  // Split the dependency strings into lists of entries
  std::vector<std::string> split_value_dependency_list =
    dealii::Utilities::split_string_list(value_dependencies);
  std::vector<std::string> split_gradient_dependency_list =
    dealii::Utilities::split_string_list(gradient_dependencies);

  // Check if either is empty and set need_value_residual and need_gradient
  // residual appropriately
  split_value_dependency_list.size() > 0 ? need_value_residual = true
                                         : need_value_residual = false;
  split_gradient_dependency_list.size() > 0 ? need_gradient_residual = true
                                            : need_gradient_residual = false;

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
      bool dependency_entry_assigned = false;

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
              need_value[dependency_variable_index] = true;
              dependency_entry_assigned             = true;

              // Check for nonlinearity
              is_nonlinear = is_nonlinear || !dependency_variable_is_explicit;
            }
          // Case if the dependency is grad(x)
          else if (dependency == gradient_variable)
            {
              need_gradient[dependency_variable_index] = true;
              dependency_entry_assigned                = true;

              // Check for nonlinearity
              is_nonlinear = is_nonlinear || !dependency_variable_is_explicit;
            }
          // Case if the dependency is hess(x)
          else if (dependency == hessian_variable)
            {
              need_hessian[dependency_variable_index] = true;
              dependency_entry_assigned               = true;

              // Check for nonlinearity
              is_nonlinear = is_nonlinear || !dependency_variable_is_explicit;
            }
          // Case if the dependency is change(x)
          else if (dependency == change_value_variable)
            {
              need_value_change[dependency_variable_index] = true;
              dependency_entry_assigned                    = true;

              Assert(variable_index != dependency_variable_index,
                     dealii::StandardExceptions::ExcMessage(
                       "PRISMS-PF Error: Dependency entry " + dependency +
                       " is not valid because the change in a variable can "
                       "only be accessed in its own governing equation."));
            }
          // Case if the dependency is grad(change(x))
          else if (dependency == change_gradient_variable)
            {
              need_gradient_change[dependency_variable_index] = true;
              dependency_entry_assigned                       = true;

              Assert(variable_index != dependency_variable_index,
                     dealii::StandardExceptions::ExcMessage(
                       "PRISMS-PF Error: Dependency entry " + dependency +
                       " is not valid because the change in a variable can "
                       "only be accessed in its own governing equation."));
            }
          // Case if the dependency is hess(change(x))
          else if (dependency == change_hessian_variable)
            {
              need_hessian_change[dependency_variable_index] = true;
              dependency_entry_assigned                      = true;

              Assert(variable_index != dependency_variable_index,
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
                                                    dependency + " is not valid."))
    }
}

void
EquationDependencyParser::pp_parse(std::vector<std::string> var_name,
                                   std::vector<std::string> pp_var_name,
                                   std::vector<std::string> sorted_dependencies_value,
                                   std::vector<std::string> sorted_dependencies_gradient)
{
  // Initialize the calculation needed flags to false
  for (unsigned int i = 0; i < var_name.size(); i++)
    {
      pp_need_value.push_back(false);
      pp_need_gradient.push_back(false);
      pp_need_hessian.push_back(false);
    }

  // Delete whitespace in the dependencies
  for (unsigned int i = 0; i < pp_var_name.size(); i++)
    {
      if (sorted_dependencies_value.size() > 0)
        {
          strip_dependency_whitespace(sorted_dependencies_value.at(i));
        }

      if (sorted_dependencies_gradient.size() > 0)
        {
          strip_dependency_whitespace(sorted_dependencies_gradient.at(i));
        }
    }

  // Now parse the dependency strings to set the flags to true where needed
  for (unsigned int i = 0; i < pp_var_name.size(); i++)
    {
      bool need_value_residual_entry, need_gradient_residual_entry;

      parseDependencyListPP(var_name,
                            sorted_dependencies_value.at(i),
                            sorted_dependencies_gradient.at(i),
                            pp_need_value,
                            pp_need_gradient,
                            pp_need_hessian,
                            need_value_residual_entry,
                            need_gradient_residual_entry);

      pp_need_value_residual.push_back(need_value_residual_entry);
      pp_need_gradient_residual.push_back(need_gradient_residual_entry);
    }
}

void
EquationDependencyParser::parseDependencyListPP(std::vector<std::string> var_name,
                                                std::string        value_dependencies,
                                                std::string        gradient_dependencies,
                                                std::vector<bool> &need_value,
                                                std::vector<bool> &need_gradient,
                                                std::vector<bool> &need_hessian,
                                                bool              &need_value_residual,
                                                bool              &need_gradient_residual)
{
  // Split the dependency strings into lists of entries
  std::vector<std::string> split_value_dependency_list =
    dealii::Utilities::split_string_list(value_dependencies);
  std::vector<std::string> split_gradient_dependency_list =
    dealii::Utilities::split_string_list(gradient_dependencies);

  // Check if either is empty and set need_value_residual and need_gradient
  // residual appropriately
  if (split_value_dependency_list.size() > 0)
    {
      need_value_residual = true;
    }
  else
    {
      need_value_residual = false;
    }

  if (split_gradient_dependency_list.size() > 0)
    {
      need_gradient_residual = true;
    }
  else
    {
      need_gradient_residual = false;
    }

  // Merge the lists of dependency entries
  /*
  std::vector<std::string> split_dependency_list = split_value_dependency_list;
  split_dependency_list.insert(split_dependency_list.end(),split_gradient_dependency_list.begin(),split_gradient_dependency_list.end());
  */

  std::vector<std::string> split_dependency_list;
  if (need_value_residual)
    {
      split_dependency_list = split_value_dependency_list;
      split_dependency_list.insert(split_dependency_list.end(),
                                   split_gradient_dependency_list.begin(),
                                   split_gradient_dependency_list.end());
    }
  else
    {
      split_dependency_list = split_gradient_dependency_list;
    }

  // Cycle through each dependency entry
  for (unsigned int dep = 0; dep < split_dependency_list.size(); dep++)
    {
      bool dependency_entry_assigned = false;

      for (unsigned int var = 0; var < var_name.size(); var++)
        {
          // Create grad() and hess() variants of the variable name
          std::string grad_var_name = {"grad()"};
          grad_var_name.insert(--grad_var_name.end(),
                               var_name.at(var).begin(),
                               var_name.at(var).end());

          std::string hess_var_name = {"hess()"};
          hess_var_name.insert(--hess_var_name.end(),
                               var_name.at(var).begin(),
                               var_name.at(var).end());

          if (split_dependency_list.at(dep) == var_name.at(var))
            {
              need_value.at(var)        = true;
              dependency_entry_assigned = true;
            }
          else if (split_dependency_list.at(dep) == grad_var_name)
            {
              need_gradient.at(var)     = true;
              dependency_entry_assigned = true;
            }
          else if (split_dependency_list.at(dep) == hess_var_name)
            {
              need_hessian.at(var)      = true;
              dependency_entry_assigned = true;
            }
        }
      if (!dependency_entry_assigned)
        {
          std::cerr << "PRISMS-PF Error: Dependency entry "
                    << split_dependency_list.at(dep) << " is not valid." << std::endl;
          abort();
        }
    }
}
