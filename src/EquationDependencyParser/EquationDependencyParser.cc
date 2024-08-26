#include "../../include/EquationDependencyParser.h"

#include <deal.II/base/utilities.h>

#include <iostream>

void
EquationDependencyParser::parse(std::vector<std::string> var_name,
                                std::vector<PDEType>     var_eq_type,
                                std::vector<std::string> sorted_dependencies_value_RHS,
                                std::vector<std::string> sorted_dependencies_gradient_RHS,
                                std::vector<std::string> sorted_dependencies_value_LHS,
                                std::vector<std::string> sorted_dependencies_gradient_LHS,
                                std::vector<bool>       &var_nonlinear)
{
  // Initialize the calculation needed flags to false
  for (unsigned int i = 0; i < var_name.size(); i++)
    {
      need_value_explicit_RHS.push_back(false);
      need_gradient_explicit_RHS.push_back(false);
      need_hessian_explicit_RHS.push_back(false);
      need_value_nonexplicit_RHS.push_back(false);
      need_gradient_nonexplicit_RHS.push_back(false);
      need_hessian_nonexplicit_RHS.push_back(false);
      need_value_nonexplicit_LHS.push_back(false);
      need_gradient_nonexplicit_LHS.push_back(false);
      need_hessian_nonexplicit_LHS.push_back(false);
      need_value_change_nonexplicit_LHS.push_back(false);
      need_gradient_change_nonexplicit_LHS.push_back(false);
      need_hessian_change_nonexplicit_LHS.push_back(false);
    }

  // Now parse the dependency strings to set the flags to true where needed
  for (unsigned int i = 0; i < var_name.size(); i++)
    {
      // First strip excess whitespace
      for (unsigned int j = 0; j < sorted_dependencies_value_RHS.at(i).length(); j++)
        {
          if (sorted_dependencies_value_RHS.at(i)[j] == ' ')
            sorted_dependencies_value_RHS.at(i).erase(j, 1);
        }
      for (unsigned int j = 0; j < sorted_dependencies_gradient_RHS.at(i).length(); j++)
        {
          if (sorted_dependencies_gradient_RHS.at(i)[j] == ' ')
            sorted_dependencies_gradient_RHS.at(i).erase(j, 1);
        }
      // Now check for each variable_eq_type
      if (var_eq_type[i] == EXPLICIT_TIME_DEPENDENT)
        {
          bool need_value_residual_entry, need_gradient_residual_entry,
            single_var_nonlinear;

          parseDependencyListRHS(var_name,
                                 var_eq_type,
                                 i,
                                 sorted_dependencies_value_RHS.at(i),
                                 sorted_dependencies_gradient_RHS.at(i),
                                 need_value_explicit_RHS,
                                 need_gradient_explicit_RHS,
                                 need_hessian_explicit_RHS,
                                 need_value_residual_entry,
                                 need_gradient_residual_entry,
                                 single_var_nonlinear);

          // std::cout << "RHS Nonlinear flag for var " << i << " :" <<
          // single_var_nonlinear << std::endl;

          var_nonlinear.push_back(single_var_nonlinear);

          need_value_residual_explicit_RHS.push_back(need_value_residual_entry);
          need_gradient_residual_explicit_RHS.push_back(need_gradient_residual_entry);

          need_value_residual_nonexplicit_RHS.push_back(false);
          need_gradient_residual_nonexplicit_RHS.push_back(false);
          need_value_residual_nonexplicit_LHS.push_back(false);
          need_gradient_residual_nonexplicit_LHS.push_back(false);
        }
      else if (var_eq_type[i] == AUXILIARY)
        {
          bool need_value_residual_entry, need_gradient_residual_entry,
            single_var_nonlinear;

          parseDependencyListRHS(var_name,
                                 var_eq_type,
                                 i,
                                 sorted_dependencies_value_RHS.at(i),
                                 sorted_dependencies_gradient_RHS.at(i),
                                 need_value_nonexplicit_RHS,
                                 need_gradient_nonexplicit_RHS,
                                 need_hessian_nonexplicit_RHS,
                                 need_value_residual_entry,
                                 need_gradient_residual_entry,
                                 single_var_nonlinear);

          var_nonlinear.push_back(single_var_nonlinear);

          // std::cout << "RHS Nonlinear flag for var " << i << " :" <<
          // single_var_nonlinear << std::endl;

          need_value_residual_explicit_RHS.push_back(false);
          need_gradient_residual_explicit_RHS.push_back(false);

          need_value_residual_nonexplicit_RHS.push_back(need_value_residual_entry);
          need_gradient_residual_nonexplicit_RHS.push_back(need_gradient_residual_entry);
          need_value_residual_nonexplicit_LHS.push_back(false);
          need_gradient_residual_nonexplicit_LHS.push_back(false);
        }
      else if (var_eq_type[i] == IMPLICIT_TIME_DEPENDENT ||
               var_eq_type[i] == TIME_INDEPENDENT)
        {
          bool need_value_residual_entry, need_gradient_residual_entry,
            single_var_nonlinear_RHS, single_var_nonlinear_LHS;

          parseDependencyListRHS(var_name,
                                 var_eq_type,
                                 i,
                                 sorted_dependencies_value_RHS.at(i),
                                 sorted_dependencies_gradient_RHS.at(i),
                                 need_value_nonexplicit_RHS,
                                 need_gradient_nonexplicit_RHS,
                                 need_hessian_nonexplicit_RHS,
                                 need_value_residual_entry,
                                 need_gradient_residual_entry,
                                 single_var_nonlinear_RHS);

          // std::cout << "RHS Nonlinear flag for var " << i << " :" <<
          // single_var_nonlinear_RHS << std::endl;

          need_value_residual_nonexplicit_RHS.push_back(need_value_residual_entry);
          need_gradient_residual_nonexplicit_RHS.push_back(need_gradient_residual_entry);

          parseDependencyListLHS(var_name,
                                 var_eq_type,
                                 i,
                                 sorted_dependencies_value_LHS.at(i),
                                 sorted_dependencies_gradient_LHS.at(i),
                                 need_value_nonexplicit_LHS,
                                 need_gradient_nonexplicit_LHS,
                                 need_hessian_nonexplicit_LHS,
                                 need_value_change_nonexplicit_LHS,
                                 need_gradient_change_nonexplicit_LHS,
                                 need_hessian_change_nonexplicit_LHS,
                                 need_value_residual_entry,
                                 need_gradient_residual_entry,
                                 single_var_nonlinear_LHS);

          // std::cout << "LHS Nonlinear flag for var " << i << " :" <<
          // single_var_nonlinear_LHS << std::endl;

          var_nonlinear.push_back(single_var_nonlinear_RHS || single_var_nonlinear_LHS);

          need_value_residual_nonexplicit_LHS.push_back(need_value_residual_entry);
          need_gradient_residual_nonexplicit_LHS.push_back(need_gradient_residual_entry);

          need_value_residual_explicit_RHS.push_back(false);
          need_gradient_residual_explicit_RHS.push_back(false);
        }
    }
}

void
EquationDependencyParser::parseDependencyListRHS(std::vector<std::string> var_name,
                                                 std::vector<PDEType>     var_eq_type,
                                                 unsigned int             var_index,
                                                 std::string        value_dependencies,
                                                 std::string        gradient_dependencies,
                                                 std::vector<bool> &need_value,
                                                 std::vector<bool> &need_gradient,
                                                 std::vector<bool> &need_hessian,
                                                 bool              &need_value_residual,
                                                 bool &need_gradient_residual,
                                                 bool &is_nonlinear)
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
  std::vector<std::string> split_dependency_list = split_value_dependency_list;
  split_dependency_list.insert(split_dependency_list.end(),
                               split_gradient_dependency_list.begin(),
                               split_gradient_dependency_list.end());

  // Cycle through each dependency entry
  // NOTE: This section is pretty confusing I think it needs refactoring or more
  // comments
  is_nonlinear = false;
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

              if ((var_eq_type[var_index] != EXPLICIT_TIME_DEPENDENT) &&
                  (var_index != var) && (var_eq_type[var] != EXPLICIT_TIME_DEPENDENT))
                {
                  is_nonlinear = true;
                }
            }
          else if (split_dependency_list.at(dep) == grad_var_name)
            {
              need_gradient.at(var)     = true;
              dependency_entry_assigned = true;
              if ((var_eq_type[var_index] != EXPLICIT_TIME_DEPENDENT) &&
                  (var_index != var) && (var_eq_type[var] != EXPLICIT_TIME_DEPENDENT))
                {
                  is_nonlinear = true;
                }
            }
          else if (split_dependency_list.at(dep) == hess_var_name)
            {
              need_hessian.at(var)      = true;
              dependency_entry_assigned = true;
              if ((var_eq_type[var_index] != EXPLICIT_TIME_DEPENDENT) &&
                  (var_index != var) && (var_eq_type[var] != EXPLICIT_TIME_DEPENDENT))
                {
                  is_nonlinear = true;
                }
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

void
EquationDependencyParser::parseDependencyListLHS(std::vector<std::string> var_name,
                                                 std::vector<PDEType>     var_eq_type,
                                                 unsigned int             var_index,
                                                 std::string        value_dependencies,
                                                 std::string        gradient_dependencies,
                                                 std::vector<bool> &need_value,
                                                 std::vector<bool> &need_gradient,
                                                 std::vector<bool> &need_hessian,
                                                 std::vector<bool> &need_value_change,
                                                 std::vector<bool> &need_gradient_change,
                                                 std::vector<bool> &need_hessian_change,
                                                 bool              &need_value_residual,
                                                 bool &need_gradient_residual,
                                                 bool &is_nonlinear)
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
  std::vector<std::string> split_dependency_list = split_value_dependency_list;
  split_dependency_list.insert(split_dependency_list.end(),
                               split_gradient_dependency_list.begin(),
                               split_gradient_dependency_list.end());

  is_nonlinear = false;
  for (unsigned int dep = 0; dep < split_dependency_list.size(); dep++)
    {
      bool dependency_entry_assigned = false;

      for (unsigned int var = 0; var < var_name.size(); var++)
        {
          std::string grad_var_name = {"grad()"};
          grad_var_name.insert(--grad_var_name.end(),
                               var_name.at(var).begin(),
                               var_name.at(var).end());

          std::string hess_var_name = {"hess()"};
          hess_var_name.insert(--hess_var_name.end(),
                               var_name.at(var).begin(),
                               var_name.at(var).end());

          std::string val_change_var_name = {"change()"};
          val_change_var_name.insert(--val_change_var_name.end(),
                                     var_name.at(var).begin(),
                                     var_name.at(var).end());

          std::string grad_change_var_name = {"grad(change())"};
          grad_change_var_name.insert(--(--grad_change_var_name.end()),
                                      var_name.at(var).begin(),
                                      var_name.at(var).end());

          std::string hess_change_var_name = {"hess(change())"};
          hess_change_var_name.insert(--(--hess_change_var_name.end()),
                                      var_name.at(var).begin(),
                                      var_name.at(var).end());

          if (split_dependency_list.at(dep) == var_name.at(var))
            {
              need_value.at(var)        = true;
              dependency_entry_assigned = true;
              if ((var_eq_type[var] != EXPLICIT_TIME_DEPENDENT))
                {
                  is_nonlinear = true;
                }
            }
          else if (split_dependency_list.at(dep) == grad_var_name)
            {
              need_gradient.at(var)     = true;
              dependency_entry_assigned = true;
              if ((var_eq_type[var] != EXPLICIT_TIME_DEPENDENT))
                {
                  is_nonlinear = true;
                }
            }
          else if (split_dependency_list.at(dep) == hess_var_name)
            {
              need_hessian.at(var)      = true;
              dependency_entry_assigned = true;
              if ((var_eq_type[var] != EXPLICIT_TIME_DEPENDENT))
                {
                  is_nonlinear = true;
                }
            }
          else if (split_dependency_list.at(dep) == val_change_var_name)
            {
              need_value_change.at(var) = true;
              dependency_entry_assigned = true;
              if (var_index != var)
                {
                  std::cerr << "PRISMS-PF Error: Dependency entry "
                            << split_dependency_list.at(dep)
                            << " is not valid because the change in a variable can "
                               "only be accessed in its own governing equation."
                            << std::endl;
                  abort();
                }
            }
          else if (split_dependency_list.at(dep) == grad_change_var_name)
            {
              need_gradient_change.at(var) = true;
              dependency_entry_assigned    = true;
              if (var_index != var)
                {
                  std::cerr << "PRISMS-PF Error: Dependency entry "
                            << split_dependency_list.at(dep)
                            << " is not valid because the change in a variable can "
                               "only be accessed in its own governing equation."
                            << std::endl;
                  abort();
                }
            }
          else if (split_dependency_list.at(dep) == hess_change_var_name)
            {
              need_hessian_change.at(var) = true;
              dependency_entry_assigned   = true;
              if (var_index != var)
                {
                  std::cerr << "PRISMS-PF Error: Dependency entry "
                            << split_dependency_list.at(dep)
                            << " is not valid because the change in a variable can "
                               "only be accessed in its own governing equation."
                            << std::endl;
                  abort();
                }
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
          for (unsigned int j = 0; j < sorted_dependencies_value.at(i).length(); j++)
            {
              if (sorted_dependencies_value.at(i)[j] == ' ')
                sorted_dependencies_value.at(i).erase(j, 1);
            }
        }

      if (sorted_dependencies_gradient.size() > 0)
        {
          for (unsigned int j = 0; j < sorted_dependencies_gradient.at(i).length(); j++)
            {
              if (sorted_dependencies_gradient.at(i)[j] == ' ')
                sorted_dependencies_gradient.at(i).erase(j, 1);
            }
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
