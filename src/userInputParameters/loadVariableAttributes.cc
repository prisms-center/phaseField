#include "../../include/userInputParameters.h"
#include "../../include/sortIndexEntryPairList.h"

template <int dim>
void userInputParameters<dim>::loadVariableAttributes(variableAttributeLoader variable_attributes){
    number_of_variables = variable_attributes.var_name_list.size();

    var_name = sortIndexEntryPairList(variable_attributes.var_name_list,number_of_variables,"var");
    var_type = sortIndexEntryPairList(variable_attributes.var_type_list,number_of_variables,SCALAR);
    var_eq_type = sortIndexEntryPairList(variable_attributes.var_eq_type_list,number_of_variables,PARABOLIC);

    // New way (v2.1) of determining the dependencies
    std::vector<std::string> sorted_dependencies_value_RHS = sortIndexEntryPairList(variable_attributes.var_eq_dependencies_value_RHS,number_of_variables,"");

    std::vector<std::string> sorted_dependencies_gradient_RHS = sortIndexEntryPairList(variable_attributes.var_eq_dependencies_gradient_RHS,number_of_variables,"");

    std::vector<std::string> sorted_dependencies_value_LHS = sortIndexEntryPairList(variable_attributes.var_eq_dependencies_value_LHS,number_of_variables,"");

    std::vector<std::string> sorted_dependencies_gradient_LHS = sortIndexEntryPairList(variable_attributes.var_eq_dependencies_gradient_LHS,number_of_variables,"");

    std::vector<bool> need_value_explicit_RHS, need_gradient_explicit_RHS, need_hessian_explicit_RHS, need_value_nonexplicit_RHS, need_gradient_nonexplicit_RHS, need_hessian_nonexplicit_RHS, need_value_nonexplicit_LHS, need_gradient_nonexplicit_LHS, need_hessian_nonexplicit_LHS, need_value_change_nonexplicit_LHS, need_gradient_change_nonexplicit_LHS, need_hessian_change_nonexplicit_LHS;

    parseDependencies(sorted_dependencies_value_RHS,
        sorted_dependencies_gradient_RHS,
        sorted_dependencies_value_LHS,
        sorted_dependencies_gradient_LHS,
        need_value_explicit_RHS,
        need_gradient_explicit_RHS,
        need_hessian_explicit_RHS,
        need_value_nonexplicit_RHS,
        need_gradient_nonexplicit_RHS,
        need_hessian_nonexplicit_RHS,
        need_value_nonexplicit_LHS,
        need_gradient_nonexplicit_LHS,
        need_hessian_nonexplicit_LHS,
        need_value_change_nonexplicit_LHS,
        need_gradient_change_nonexplicit_LHS,
        need_hessian_change_nonexplicit_LHS);

    // Sort the variable attributes and load them into individual vectors
    std::vector<bool> nucleating_variable = sortIndexEntryPairList(variable_attributes.nucleating_variable_list,number_of_variables,false);
    std::vector<bool> need_value_nucleation = sortIndexEntryPairList(variable_attributes.need_value_list_nucleation,number_of_variables,false);

    std::vector<bool> need_value = sortIndexEntryPairList(variable_attributes.need_value_list,number_of_variables,true);
	std::vector<bool> need_gradient = sortIndexEntryPairList(variable_attributes.need_gradient_list,number_of_variables,true);
	std::vector<bool> need_hessian = sortIndexEntryPairList(variable_attributes.need_hessian_list,number_of_variables,false);
	std::vector<bool> value_residual = sortIndexEntryPairList(variable_attributes.need_value_residual_list,number_of_variables,false);
	std::vector<bool> gradient_residual = sortIndexEntryPairList(variable_attributes.need_gradient_residual_list,number_of_variables,false);

    std::vector<bool> need_value_LHS = sortIndexEntryPairList(variable_attributes.need_value_list_LHS,number_of_variables,false);
	std::vector<bool> need_gradient_LHS = sortIndexEntryPairList(variable_attributes.need_gradient_list_LHS,number_of_variables,false);
	std::vector<bool> need_hessian_LHS = sortIndexEntryPairList(variable_attributes.need_hessian_list_LHS,number_of_variables,false);
	std::vector<bool> value_residual_LHS = sortIndexEntryPairList(variable_attributes.need_value_residual_list_LHS,number_of_variables,false);
	std::vector<bool> gradient_residual_LHS = sortIndexEntryPairList(variable_attributes.need_gradient_residual_list_LHS,number_of_variables,false);

    std::vector<bool> need_value_change_LHS = sortIndexEntryPairList(variable_attributes.need_value_change_list_LHS,number_of_variables,false);
	std::vector<bool> need_gradient_change_LHS = sortIndexEntryPairList(variable_attributes.need_gradient_change_list_LHS,number_of_variables,false);
	std::vector<bool> need_hessian_change_LHS = sortIndexEntryPairList(variable_attributes.need_hessian_change_list_LHS,number_of_variables,false);

    std::vector<bool> need_value_pp = sortIndexEntryPairList(variable_attributes.need_value_list_PP,number_of_variables,true);
	std::vector<bool> need_gradient_pp = sortIndexEntryPairList(variable_attributes.need_gradient_list_PP,number_of_variables,true);
	std::vector<bool> need_hessian_pp = sortIndexEntryPairList(variable_attributes.need_hessian_list_PP,number_of_variables,true);

    // Sort the post-processing variable attributes
    pp_number_of_variables = variable_attributes.var_name_list_PP.size();
    pp_var_name = sortIndexEntryPairList(variable_attributes.var_name_list_PP,pp_number_of_variables,"var");
    pp_var_type = sortIndexEntryPairList(variable_attributes.var_type_list_PP,pp_number_of_variables,SCALAR);
    std::vector<bool> pp_value_residual = sortIndexEntryPairList(variable_attributes.need_value_residual_list_PP,pp_number_of_variables,false);
	std::vector<bool> pp_gradient_residual = sortIndexEntryPairList(variable_attributes.need_gradient_residual_list_PP,pp_number_of_variables,false);
    pp_calc_integral = sortIndexEntryPairList(variable_attributes.output_integral_list,pp_number_of_variables,false);


    // Load some nucleation parameters
    for (unsigned int i=0; i<number_of_variables; i++){
        if (nucleating_variable.at(i)==true){
            nucleating_variable_indices.push_back(i);
        }
        if (need_value_nucleation.at(i) || nucleating_variable.at(i)){
            nucleation_need_value.push_back(i);
        }
    }

    if (nucleating_variable_indices.size() > 0){
        nucleation_occurs = true;
    }
    else {
        nucleation_occurs = false;
    }


    // Load these attributes into the varInfoList objects
	varInfoListRHS.reserve(number_of_variables);
	unsigned int scalar_var_index = 0;
	unsigned int vector_var_index = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		variable_info varInfo;

        varInfo.need_value = need_value[i];
        varInfo.need_gradient = need_gradient[i];
        varInfo.need_hessian = need_hessian[i];
        varInfo.value_residual = value_residual[i];
        varInfo.gradient_residual = gradient_residual[i];

        varInfo.global_var_index = i;

		if (need_value[i] or need_gradient[i] or need_hessian[i]){
            varInfo.var_needed = true;
		}
        else {
            varInfo.var_needed = false;
        }

        if (var_type[i] == SCALAR){
            varInfo.is_scalar = true;
            if (varInfo.var_needed){
                varInfo.scalar_or_vector_index = scalar_var_index;
                scalar_var_index++;
            }
        }
        else {
            varInfo.is_scalar = false;
            if (varInfo.var_needed){
                varInfo.scalar_or_vector_index = vector_var_index;
                vector_var_index++;
            }
        }

        varInfoListRHS.push_back(varInfo);
	}

	// Load variable information for calculating the LHS
	num_var_LHS = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		if (need_value_LHS[i] or need_gradient_LHS[i] or need_hessian_LHS[i]){
			num_var_LHS++;
		}
	}

	varInfoListLHS.reserve(num_var_LHS);
	scalar_var_index = 0;
	vector_var_index = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		variable_info varInfo;

        varInfo.need_value = need_value_LHS[i];
        varInfo.need_gradient = need_gradient_LHS[i];
        varInfo.need_hessian = need_hessian_LHS[i];
        varInfo.value_residual = value_residual_LHS[i];
        varInfo.gradient_residual = gradient_residual_LHS[i];

        varInfo.global_var_index = i;


		if (need_value_LHS[i] or need_gradient_LHS[i] or need_hessian_LHS[i]){
            varInfo.var_needed = true;
		}
        else {
            varInfo.var_needed = false;
		}

        if (var_type[i] == SCALAR){
            varInfo.is_scalar = true;
            if (varInfo.var_needed){
                varInfo.scalar_or_vector_index = scalar_var_index;
                scalar_var_index++;
            }
        }
        else {
            varInfo.is_scalar = false;
            if (varInfo.var_needed){
                varInfo.scalar_or_vector_index = vector_var_index;
                vector_var_index++;
            }
        }

        varInfoListLHS.push_back(varInfo);
	}

    varChangeInfoListLHS.reserve(num_var_LHS);
	scalar_var_index = 0;
	vector_var_index = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		variable_info varInfo;

        varInfo.need_value = need_value_change_LHS[i];
        varInfo.need_gradient = need_gradient_change_LHS[i];
        varInfo.need_hessian = need_hessian_change_LHS[i];

        // FOR NOW, TAKING THESE FROM THE VARIABLE ITSELF!!
        varInfo.value_residual = value_residual_LHS[i];
        varInfo.gradient_residual = gradient_residual_LHS[i];

        varInfo.global_var_index = i;

		if (varInfo.need_value or varInfo.need_gradient or varInfo.need_hessian){
            varInfo.var_needed = true;
		}
        else {
            varInfo.var_needed = false;
		}

        if (var_type[i] == SCALAR){
            varInfo.is_scalar = true;
            if (varInfo.var_needed){
                varInfo.scalar_or_vector_index = scalar_var_index;
                scalar_var_index++;
            }
        }
        else {
            varInfo.is_scalar = false;
            if (varInfo.var_needed){
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
	for (unsigned int i=0; i<number_of_variables; i++){
		variable_info varInfo;

        varInfo.need_value = need_value_pp[i];
        varInfo.need_gradient = need_gradient_pp[i];
        varInfo.need_hessian = need_hessian_pp[i];

        varInfo.global_var_index = i;

		if (need_value_pp[i] or need_gradient_pp[i] or need_hessian_pp[i]){
            varInfo.var_needed = true;
		}
        else {
            varInfo.var_needed = false;
        }

        if (var_type[i] == SCALAR){
            varInfo.is_scalar = true;
            if (varInfo.var_needed){
                varInfo.scalar_or_vector_index = scalar_var_index;
                scalar_var_index++;
            }
        }
        else {
            varInfo.is_scalar = false;
            if (varInfo.var_needed){
                varInfo.scalar_or_vector_index = vector_var_index;
                vector_var_index++;
            }
        }

        pp_baseVarInfoList.push_back(varInfo);
	}

    // Now load the information for the post-processing variables
    // Parameters for postprocessing
    if (pp_number_of_variables > 0){
        postProcessingRequired = true;
    }
    else {
        postProcessingRequired = false;
    }

    num_integrated_fields = 0;
    for (unsigned int i=0; i<pp_number_of_variables; i++){
        if (pp_calc_integral[i]){
            num_integrated_fields++;
            integrated_field_indices.push_back(i);
        }
    }

    // The info list for the postprocessing field variables
	pp_varInfoList.reserve(pp_number_of_variables);
	scalar_var_index = 0;
	vector_var_index = 0;
	for (unsigned int i=0; i<pp_number_of_variables; i++){
		variable_info varInfo;
        varInfo.var_needed = true;

        varInfo.value_residual = pp_value_residual[i];
        varInfo.gradient_residual = pp_gradient_residual[i];

		varInfo.global_var_index = i;
		if (pp_var_type[i] == SCALAR){
			varInfo.is_scalar = true;
			varInfo.scalar_or_vector_index = scalar_var_index;
			scalar_var_index++;
		}
		else {
			varInfo.is_scalar = false;
			varInfo.scalar_or_vector_index = vector_var_index;
			vector_var_index++;
		}
		pp_varInfoList.push_back(varInfo);
    }

}

// Template instantiations
#include "../../include/userInputParameters_template_instantiations.h"
