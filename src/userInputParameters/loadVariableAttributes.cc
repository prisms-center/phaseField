#include "../../include/EquationDependencyParser.h"
#include "../../include/userInputParameters.h"

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

  nucleating_variable_indices.size() > 0 ? nucleation_occurs = true
                                         : nucleation_occurs = false;

  // Load these attributes into the varInfoList objects

  // Load variable information for calculating the RHS for explicit equations
  num_var_explicit_RHS = 0;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      if (!(variable_attributes.equation_dependency_parser.eval_flags_explicit_RHS[i] &
            dealii::EvaluationFlags::nothing))
        {
          num_var_explicit_RHS++;
        }
    }
  varInfoListExplicitRHS.reserve(num_var_explicit_RHS);
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      variable_info varInfo;

      varInfo.evaluation_flags =
        variable_attributes.equation_dependency_parser.eval_flags_explicit_RHS[i];

      varInfo.residual_flags = variable_attributes.equation_dependency_parser
                                 .eval_flags_residual_explicit_RHS[i];

      varInfo.global_var_index = i;

      varInfo.var_needed = !(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = var_type[i] == SCALAR;

      varInfoListExplicitRHS.push_back(varInfo);
    }

  // Load variable information for calculating the RHS for nonexplicit equations
  num_var_nonexplicit_RHS = 0;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      if (!(variable_attributes.equation_dependency_parser.eval_flags_nonexplicit_RHS[i] &
            dealii::EvaluationFlags::nothing))
        {
          num_var_nonexplicit_RHS++;
        }
    }
  varInfoListNonexplicitRHS.reserve(num_var_nonexplicit_RHS);
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      variable_info varInfo;

      varInfo.evaluation_flags =
        variable_attributes.equation_dependency_parser.eval_flags_nonexplicit_RHS[i];

      varInfo.residual_flags = variable_attributes.equation_dependency_parser
                                 .eval_flags_residual_nonexplicit_RHS[i];

      varInfo.global_var_index = i;

      varInfo.var_needed = !(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = var_type[i] == SCALAR;

      varInfoListNonexplicitRHS.push_back(varInfo);
    }

  // Load variable information for calculating the LHS
  num_var_LHS = 0;
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      if (!(variable_attributes.equation_dependency_parser.eval_flags_nonexplicit_LHS[i] &
            dealii::EvaluationFlags::nothing))
        {
          num_var_LHS++;
        }
    }

  varInfoListLHS.reserve(num_var_LHS);
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      variable_info varInfo;

      varInfo.evaluation_flags =
        variable_attributes.equation_dependency_parser.eval_flags_nonexplicit_LHS[i];

      varInfo.residual_flags = variable_attributes.equation_dependency_parser
                                 .eval_flags_residual_nonexplicit_LHS[i];

      varInfo.global_var_index = i;

      varInfo.var_needed = !(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = var_type[i] == SCALAR;

      varInfoListLHS.push_back(varInfo);
    }

  varChangeInfoListLHS.reserve(num_var_LHS);
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      variable_info varInfo;

      varInfo.evaluation_flags = variable_attributes.equation_dependency_parser
                                   .eval_flags_change_nonexplicit_LHS[i];

      // FOR NOW, TAKING THESE FROM THE VARIABLE ITSELF!!
      varInfo.residual_flags = variable_attributes.equation_dependency_parser
                                 .eval_flags_residual_nonexplicit_LHS[i];

      varInfo.global_var_index = i;

      varInfo.var_needed = !(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = var_type[i] == SCALAR;

      varChangeInfoListLHS.push_back(varInfo);
    }

  // Load variable information for postprocessing
  // First, the info list for the base field variables
  pp_baseVarInfoList.reserve(number_of_variables);
  for (unsigned int i = 0; i < number_of_variables; i++)
    {
      variable_info varInfo;

      varInfo.evaluation_flags =
        variable_attributes.equation_dependency_parser.eval_flags_postprocess[i];

      varInfo.global_var_index = i;

      varInfo.var_needed = !(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = var_type[i] == SCALAR;

      pp_baseVarInfoList.push_back(varInfo);
    }

  // Now load the information for the post-processing variables
  // Parameters for postprocessing
  pp_number_of_variables > 0 ? postProcessingRequired = true
                             : postProcessingRequired = false;

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
  for (unsigned int i = 0; i < pp_number_of_variables; i++)
    {
      variable_info varInfo;
      varInfo.var_needed = true;

      varInfo.residual_flags =
        variable_attributes.equation_dependency_parser.eval_flags_residual_postprocess[i];

      varInfo.global_var_index = i;

      varInfo.is_scalar = pp_var_type[i] == SCALAR;

      pp_varInfoList.push_back(varInfo);
    }
}

// Template instantiations
#include "../../include/userInputParameters_template_instantiations.h"
