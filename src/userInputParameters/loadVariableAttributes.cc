#include "../../include/userInputParameters.h"

template <int dim>
void
userInputParameters<dim>::loadVariableAttributes(
  const variableAttributeLoader &variable_attributes)
{
  number_of_variables    = variable_attributes.attributes.size();
  pp_number_of_variables = variable_attributes.pp_attributes.size();
  // Load some nucleation parameters
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      if (variable.nucleating_variable)
        {
          nucleating_variable_indices.push_back(index);
        }
      if (variable.need_value_nucleation || variable.nucleating_variable)
        {
          nucleation_need_value.push_back(index);
        }
    }

  nucleation_occurs = !nucleating_variable_indices.empty();

  // Load variable information for calculating the RHS for explicit equations
  num_var_explicit_RHS = 0;
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      if (!(variable.eval_flags_explicit_RHS & dealii::EvaluationFlags::nothing))
        {
          num_var_explicit_RHS++;
        }
    }
  varInfoListExplicitRHS.reserve(num_var_explicit_RHS);
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      variable_info varInfo;

      varInfo.evaluation_flags = variable.eval_flags_explicit_RHS;

      varInfo.residual_flags = variable.eval_flags_residual_explicit_RHS;

      varInfo.global_var_index = index;

      varInfo.var_needed = !(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = variable.var_type == SCALAR;

      varInfoListExplicitRHS.push_back(varInfo);
    }

  // Load variable information for calculating the RHS for nonexplicit equations
  num_var_nonexplicit_RHS = 0;
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      if (!(variable.eval_flags_nonexplicit_RHS & dealii::EvaluationFlags::nothing))
        {
          num_var_nonexplicit_RHS++;
        }
    }
  varInfoListNonexplicitRHS.reserve(num_var_nonexplicit_RHS);
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      variable_info varInfo;

      varInfo.evaluation_flags = variable.eval_flags_nonexplicit_RHS;

      varInfo.residual_flags = variable.eval_flags_residual_nonexplicit_RHS;

      varInfo.global_var_index = index;

      varInfo.var_needed = !(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = variable.var_type == SCALAR;

      varInfoListNonexplicitRHS.push_back(varInfo);
    }

  // Load variable information for calculating the LHS
  num_var_LHS = 0;
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      if (!(variable.eval_flags_nonexplicit_LHS & dealii::EvaluationFlags::nothing))
        {
          num_var_LHS++;
        }
    }

  varInfoListLHS.reserve(num_var_LHS);
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      variable_info varInfo;

      varInfo.evaluation_flags = variable.eval_flags_nonexplicit_LHS;

      varInfo.residual_flags = variable.eval_flags_residual_nonexplicit_LHS;

      varInfo.global_var_index = index;

      varInfo.var_needed = !(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = variable.var_type == SCALAR;

      varInfoListLHS.push_back(varInfo);
    }

  varChangeInfoListLHS.reserve(num_var_LHS);
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      variable_info varInfo;

      varInfo.evaluation_flags = variable.eval_flags_change_nonexplicit_LHS;

      // FOR NOW, TAKING THESE FROM THE VARIABLE ITSELF!!
      varInfo.residual_flags = variable.eval_flags_residual_nonexplicit_LHS;

      varInfo.global_var_index = index;

      varInfo.var_needed = !(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = variable.var_type == SCALAR;

      varChangeInfoListLHS.push_back(varInfo);
    }

  // Load variable information for postprocessing
  // First, the info list for the base field variables
  pp_baseVarInfoList.reserve(number_of_variables);
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      variable_info varInfo;

      varInfo.evaluation_flags = variable.eval_flags_postprocess;

      varInfo.global_var_index = index;

      varInfo.var_needed = !(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = variable.var_type == SCALAR;

      pp_baseVarInfoList.push_back(varInfo);
    }

  // Now load the information for the post-processing variables
  // Parameters for postprocessing

  postProcessingRequired = pp_number_of_variables > 0;

  num_integrated_fields = 0;
  for (const auto &[pp_index, pp_variable] : variable_attributes.pp_attributes)
    {
      if (pp_variable.calc_integral)
        {
          num_integrated_fields++;
          integrated_field_indices.push_back(pp_index);
        }
    }

  // The info list for the postprocessing field variables
  pp_varInfoList.reserve(pp_number_of_variables);
  for (const auto &[pp_index, pp_variable] : variable_attributes.pp_attributes)
    {
      variable_info varInfo;

      varInfo.var_needed = true;

      varInfo.residual_flags = pp_variable.eval_flags_residual_postprocess;

      varInfo.global_var_index = pp_index;

      varInfo.is_scalar = pp_variable.var_type == SCALAR;

      pp_varInfoList.push_back(varInfo);
    }
}

// Template instantiations
#include "../../include/userInputParameters_template_instantiations.h"
