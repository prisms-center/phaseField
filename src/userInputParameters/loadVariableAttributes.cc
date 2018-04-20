#include "../../include/userInputParameters.h"
#include "../../include/sortIndexEntryPairList.h"
#include "../../include/EquationDependencyParser.h"

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

    EquationDependencyParser equation_dependency_parser(
        var_name,
        var_eq_type,
        sorted_dependencies_value_RHS,
        sorted_dependencies_gradient_RHS,
        sorted_dependencies_value_LHS,
        sorted_dependencies_gradient_LHS, var_nonlinear);

    // NOTE: THIS NEEDS TO GO INSIDE EquationDependencyParser
    // Determine which equations are nonlinear
    // Equations are nonlinear if the LHS has non-change values of the equation's variable or if it is
    // ellptic and uses variables from other elliptic equations
    for (unsigned int i=0; i<var_name.size(); i++){
        var_nonlinear.push_back(false);
        if (var_eq_type[i] == ELLIPTIC){
            if (equation_dependency_parser.need_value_nonexplicit_LHS[i] or equation_dependency_parser.need_gradient_nonexplicit_LHS[i] or equation_dependency_parser.need_hessian_nonexplicit_LHS[i]){
                var_nonlinear.at(i) = true;
            }
            else {
                for (unsigned int j=0; j<var_name.size(); j++){
                    if (var_eq_type[j] == ELLIPTIC and j != i ){
                        if (equation_dependency_parser.need_value_nonexplicit_RHS[i] or equation_dependency_parser.need_gradient_nonexplicit_RHS[i] or equation_dependency_parser.need_hessian_nonexplicit_RHS[i] or equation_dependency_parser.need_value_nonexplicit_LHS[i] or equation_dependency_parser.need_gradient_nonexplicit_LHS[i] or equation_dependency_parser.need_hessian_nonexplicit_LHS[i])
                    }
                }


                var_nonlinear.at(i) = true;
            }
        }
    }


    // Sort the variable attributes and load them into individual vectors
    std::vector<bool> nucleating_variable = sortIndexEntryPairList(variable_attributes.nucleating_variable_list,number_of_variables,false);
    std::vector<bool> need_value_nucleation = sortIndexEntryPairList(variable_attributes.need_value_list_nucleation,number_of_variables,false);

    /*
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
    */

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

    // Load variable information for calculating the RHS for explicit equations
    num_var_explicit_RHS = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		if (equation_dependency_parser.need_value_explicit_RHS[i] or equation_dependency_parser.need_gradient_explicit_RHS[i] or equation_dependency_parser.need_hessian_explicit_RHS[i]){
			num_var_explicit_RHS++;
		}
	}
	varInfoListExplicitRHS.reserve(num_var_explicit_RHS);
	unsigned int scalar_var_index = 0;
	unsigned int vector_var_index = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		variable_info varInfo;

        varInfo.need_value = equation_dependency_parser.need_value_explicit_RHS[i];
        varInfo.need_gradient = equation_dependency_parser.need_gradient_explicit_RHS[i];
        varInfo.need_hessian = equation_dependency_parser.need_hessian_explicit_RHS[i];
        varInfo.value_residual = equation_dependency_parser.need_value_residual_explicit_RHS[i];
        varInfo.gradient_residual = equation_dependency_parser.need_gradient_residual_explicit_RHS[i];

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

        varInfoListExplicitRHS.push_back(varInfo);
	}

    // Load variable information for calculating the RHS for nonexplicit equations
    num_var_nonexplicit_RHS = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		if (equation_dependency_parser.need_value_nonexplicit_RHS[i] or equation_dependency_parser.need_gradient_nonexplicit_RHS[i] or equation_dependency_parser.need_hessian_nonexplicit_RHS[i]){
			num_var_nonexplicit_RHS++;
		}
	}
	varInfoListNonexplicitRHS.reserve(num_var_nonexplicit_RHS);
	scalar_var_index = 0;
	vector_var_index = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		variable_info varInfo;

        varInfo.need_value = equation_dependency_parser.need_value_nonexplicit_RHS[i];
        varInfo.need_gradient = equation_dependency_parser.need_gradient_nonexplicit_RHS[i];
        varInfo.need_hessian = equation_dependency_parser.need_hessian_nonexplicit_RHS[i];
        varInfo.value_residual = equation_dependency_parser.need_value_residual_nonexplicit_RHS[i];
        varInfo.gradient_residual = equation_dependency_parser.need_gradient_residual_nonexplicit_RHS[i];

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

        varInfoListNonexplicitRHS.push_back(varInfo);
	}

	// Load variable information for calculating the LHS
	num_var_LHS = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		if (equation_dependency_parser.need_value_nonexplicit_LHS[i] or equation_dependency_parser.need_gradient_nonexplicit_LHS[i] or equation_dependency_parser.need_hessian_nonexplicit_LHS[i]){
			num_var_LHS++;
		}
	}

	varInfoListLHS.reserve(num_var_LHS);
	scalar_var_index = 0;
	vector_var_index = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		variable_info varInfo;

        varInfo.need_value = equation_dependency_parser.need_value_nonexplicit_LHS[i];
        varInfo.need_gradient = equation_dependency_parser.need_gradient_nonexplicit_LHS[i];
        varInfo.need_hessian = equation_dependency_parser.need_hessian_nonexplicit_LHS[i];
        varInfo.value_residual = equation_dependency_parser.need_value_residual_nonexplicit_LHS[i];
        varInfo.gradient_residual = equation_dependency_parser.need_gradient_residual_nonexplicit_LHS[i];

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

        varInfoListLHS.push_back(varInfo);
	}

    varChangeInfoListLHS.reserve(num_var_LHS);
	scalar_var_index = 0;
	vector_var_index = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		variable_info varInfo;

        varInfo.need_value = equation_dependency_parser.need_value_change_nonexplicit_LHS[i];
        varInfo.need_gradient = equation_dependency_parser.need_gradient_change_nonexplicit_LHS[i];
        varInfo.need_hessian = equation_dependency_parser.need_hessian_change_nonexplicit_LHS[i];

        // FOR NOW, TAKING THESE FROM THE VARIABLE ITSELF!!
        varInfo.value_residual = equation_dependency_parser.need_value_residual_nonexplicit_LHS[i];
        varInfo.gradient_residual = equation_dependency_parser.need_gradient_residual_nonexplicit_LHS[i];

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
