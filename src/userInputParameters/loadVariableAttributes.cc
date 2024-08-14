#include "../../include/userInputParameters.h"
// #include "../../include/sortIndexEntryPairList.h"
#include "../../include/EquationDependencyParser.h"

template <int dim>
void
userInputParameters<dim>::loadVariableAttributes(
  variableAttributeLoader variable_attributes)
{
  // Pull some variable values from variable_attributes
  number_of_variables = variable_attributes.number_of_variables;
  var_name            = variable_attributes.var_name;
  var_type            = variable_attributes.var_type;
  var_eq_type         = variable_attributes.var_eq_type;

  var_nonlinear = variable_attributes.var_nonlinear;

  pp_calc_integral = variable_attributes.pp_calc_integral;

  pp_number_of_variables = variable_attributes.pp_number_of_variables;
  pp_var_name            = variable_attributes.pp_var_name;
  pp_var_type            = variable_attributes.pp_var_type;

  // Load some nucleation parameters
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      if (variable_attributes.nucleating_variable.at(i) == true)
        {
          nucleating_variable_indices.push_back(i);
        }
      if (variable_attributes.need_value_nucleation.at(i) ||
          variable_attributes.nucleating_variable.at(i))
        {
          nucleation_need_value.push_back(i);
        }
    }

  if (nucleating_variable_indices.size() > 0)
    {
      nucleation_occurs = true;
    }
  else
    {
      nucleation_occurs = false;
    }

  // Load these attributes into the varInfoList objects

  // Load variable information for calculating the RHS for explicit equations
  num_var_explicit_RHS = 0;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      if (variable_attributes.equation_dependency_parser.need_value_explicit_RHS[i] or
          variable_attributes.equation_dependency_parser.need_gradient_explicit_RHS[i] or
          variable_attributes.equation_dependency_parser.need_hessian_explicit_RHS[i])
        {
          num_var_explicit_RHS++;
        }
    }
  varInfoListExplicitRHS.reserve(num_var_explicit_RHS);
  unsigned int scalar_var_index = 0;
  unsigned int vector_var_index = 0;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      variable_info varInfo;

      varInfo.need_value =
        variable_attributes.equation_dependency_parser.need_value_explicit_RHS[i];
      varInfo.need_gradient =
        variable_attributes.equation_dependency_parser.need_gradient_explicit_RHS[i];
      varInfo.need_hessian =
        variable_attributes.equation_dependency_parser.need_hessian_explicit_RHS[i];
      varInfo.value_residual = variable_attributes.equation_dependency_parser
                                 .need_value_residual_explicit_RHS[i];
      varInfo.gradient_residual = variable_attributes.equation_dependency_parser
                                    .need_gradient_residual_explicit_RHS[i];

      varInfo.global_var_index = i;

      if (varInfo.need_value or varInfo.need_gradient or varInfo.need_hessian)
        {
          varInfo.var_needed = true;
        }
      else
        {
          varInfo.var_needed = false;
        }

      if (var_type[i] == SCALAR)
        {
          varInfo.is_scalar = true;
          if (varInfo.var_needed)
            {
              varInfo.scalar_or_vector_index = scalar_var_index;
              scalar_var_index++;
            }
        }
      else
        {
          varInfo.is_scalar = false;
          if (varInfo.var_needed)
            {
              varInfo.scalar_or_vector_index = vector_var_index;
              vector_var_index++;
            }
        }

      varInfoListExplicitRHS.push_back(varInfo);
    }

  // Load variable information for calculating the RHS for nonexplicit equations
  num_var_nonexplicit_RHS = 0;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      if (variable_attributes.equation_dependency_parser.need_value_nonexplicit_RHS[i] or
          variable_attributes.equation_dependency_parser
            .need_gradient_nonexplicit_RHS[i] or
          variable_attributes.equation_dependency_parser.need_hessian_nonexplicit_RHS[i])
        {
          num_var_nonexplicit_RHS++;
        }
    }
  varInfoListNonexplicitRHS.reserve(num_var_nonexplicit_RHS);
  scalar_var_index = 0;
  vector_var_index = 0;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      variable_info varInfo;

      varInfo.need_value =
        variable_attributes.equation_dependency_parser.need_value_nonexplicit_RHS[i];
      varInfo.need_gradient =
        variable_attributes.equation_dependency_parser.need_gradient_nonexplicit_RHS[i];
      varInfo.need_hessian =
        variable_attributes.equation_dependency_parser.need_hessian_nonexplicit_RHS[i];
      varInfo.value_residual = variable_attributes.equation_dependency_parser
                                 .need_value_residual_nonexplicit_RHS[i];
      varInfo.gradient_residual = variable_attributes.equation_dependency_parser
                                    .need_gradient_residual_nonexplicit_RHS[i];

      varInfo.global_var_index = i;

      if (varInfo.need_value or varInfo.need_gradient or varInfo.need_hessian)
        {
          varInfo.var_needed = true;
        }
      else
        {
          varInfo.var_needed = false;
        }

      if (var_type[i] == SCALAR)
        {
          varInfo.is_scalar = true;
          if (varInfo.var_needed)
            {
              varInfo.scalar_or_vector_index = scalar_var_index;
              scalar_var_index++;
            }
        }
      else
        {
          varInfo.is_scalar = false;
          if (varInfo.var_needed)
            {
              varInfo.scalar_or_vector_index = vector_var_index;
              vector_var_index++;
            }
        }

      varInfoListNonexplicitRHS.push_back(varInfo);
    }

  // Load variable information for calculating the LHS
  num_var_LHS = 0;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      if (variable_attributes.equation_dependency_parser.need_value_nonexplicit_LHS[i] or
          variable_attributes.equation_dependency_parser
            .need_gradient_nonexplicit_LHS[i] or
          variable_attributes.equation_dependency_parser.need_hessian_nonexplicit_LHS[i])
        {
          num_var_LHS++;
        }
    }

  varInfoListLHS.reserve(num_var_LHS);
  scalar_var_index = 0;
  vector_var_index = 0;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      variable_info varInfo;

      varInfo.need_value =
        variable_attributes.equation_dependency_parser.need_value_nonexplicit_LHS[i];
      varInfo.need_gradient =
        variable_attributes.equation_dependency_parser.need_gradient_nonexplicit_LHS[i];
      varInfo.need_hessian =
        variable_attributes.equation_dependency_parser.need_hessian_nonexplicit_LHS[i];
      varInfo.value_residual = variable_attributes.equation_dependency_parser
                                 .need_value_residual_nonexplicit_LHS[i];
      varInfo.gradient_residual = variable_attributes.equation_dependency_parser
                                    .need_gradient_residual_nonexplicit_LHS[i];

      varInfo.global_var_index = i;

      if (varInfo.need_value or varInfo.need_gradient or varInfo.need_hessian)
        {
          varInfo.var_needed = true;
        }
      else
        {
          varInfo.var_needed = false;
        }

      if (var_type[i] == SCALAR)
        {
          varInfo.is_scalar = true;
          if (varInfo.var_needed)
            {
              varInfo.scalar_or_vector_index = scalar_var_index;
              scalar_var_index++;
            }
        }
      else
        {
          varInfo.is_scalar = false;
          if (varInfo.var_needed)
            {
              varInfo.scalar_or_vector_index = vector_var_index;
              vector_var_index++;
            }
        }

      varInfoListLHS.push_back(varInfo);
    }

  varChangeInfoListLHS.reserve(num_var_LHS);
  scalar_var_index = 0;
  vector_var_index = 0;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      variable_info varInfo;

      varInfo.need_value = variable_attributes.equation_dependency_parser
                             .need_value_change_nonexplicit_LHS[i];
      varInfo.need_gradient = variable_attributes.equation_dependency_parser
                                .need_gradient_change_nonexplicit_LHS[i];
      varInfo.need_hessian = variable_attributes.equation_dependency_parser
                               .need_hessian_change_nonexplicit_LHS[i];

      // FOR NOW, TAKING THESE FROM THE VARIABLE ITSELF!!
      varInfo.value_residual = variable_attributes.equation_dependency_parser
                                 .need_value_residual_nonexplicit_LHS[i];
      varInfo.gradient_residual = variable_attributes.equation_dependency_parser
                                    .need_gradient_residual_nonexplicit_LHS[i];

      varInfo.global_var_index = i;

      if (varInfo.need_value or varInfo.need_gradient or varInfo.need_hessian)
        {
          varInfo.var_needed = true;
        }
      else
        {
          varInfo.var_needed = false;
        }

      if (var_type[i] == SCALAR)
        {
          varInfo.is_scalar = true;
          if (varInfo.var_needed)
            {
              varInfo.scalar_or_vector_index = scalar_var_index;
              scalar_var_index++;
            }
        }
      else
        {
          varInfo.is_scalar = false;
          if (varInfo.var_needed)
            {
              varInfo.scalar_or_vector_index = vector_var_index;
              vector_var_index++;
            }
        }

      varChangeInfoListLHS.push_back(varInfo);
    }

  // Load variable information for postprocessing
  // First, the info list for the base field variables
  pp_baseVarInfoList.reserve(number_of_variables);
  scalar_var_index = 0;
  vector_var_index = 0;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      variable_info varInfo;

      varInfo.need_value =
        variable_attributes.equation_dependency_parser.pp_need_value[i];
      varInfo.need_gradient =
        variable_attributes.equation_dependency_parser.pp_need_gradient[i];
      varInfo.need_hessian =
        variable_attributes.equation_dependency_parser.pp_need_hessian[i];

      varInfo.global_var_index = i;

      if (variable_attributes.equation_dependency_parser.pp_need_value[i] or
          variable_attributes.equation_dependency_parser.pp_need_gradient[i] or
          variable_attributes.equation_dependency_parser.pp_need_hessian[i])
        {
          varInfo.var_needed = true;
        }
      else
        {
          varInfo.var_needed = false;
        }

      if (var_type[i] == SCALAR)
        {
          varInfo.is_scalar = true;
          if (varInfo.var_needed)
            {
              varInfo.scalar_or_vector_index = scalar_var_index;
              scalar_var_index++;
            }
        }
      else
        {
          varInfo.is_scalar = false;
          if (varInfo.var_needed)
            {
              varInfo.scalar_or_vector_index = vector_var_index;
              vector_var_index++;
            }
        }

      pp_baseVarInfoList.push_back(varInfo);
    }

  // Now load the information for the post-processing variables
  // Parameters for postprocessing
  if (pp_number_of_variables > 0)
    {
      postProcessingRequired = true;
    }
  else
    {
      postProcessingRequired = false;
    }

  num_integrated_fields = 0;
  for (unsigned int i = 0; i < pp_number_of_variables; i++)
    {
      if (pp_calc_integral[i])
        {
          num_integrated_fields++;
          integrated_field_indices.push_back(i);
        }
    }

  // The info list for the postprocessing field variables
  pp_varInfoList.reserve(pp_number_of_variables);
  scalar_var_index = 0;
  vector_var_index = 0;
  for (unsigned int i = 0; i < pp_number_of_variables; i++)
    {
      variable_info varInfo;
      varInfo.var_needed = true;

      varInfo.value_residual =
        variable_attributes.equation_dependency_parser.pp_need_value_residual[i];
      varInfo.gradient_residual =
        variable_attributes.equation_dependency_parser.pp_need_gradient_residual[i];

      varInfo.global_var_index = i;
      if (pp_var_type[i] == SCALAR)
        {
          varInfo.is_scalar              = true;
          varInfo.scalar_or_vector_index = scalar_var_index;
          scalar_var_index++;
        }
      else
        {
          varInfo.is_scalar              = false;
          varInfo.scalar_or_vector_index = vector_var_index;
          vector_var_index++;
        }
      pp_varInfoList.push_back(varInfo);
    }
}

// Template instantiations
#include "../../include/userInputParameters_template_instantiations.h"
