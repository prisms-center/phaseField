#include "../../include/dealIIheaders.h"
#include "../../include/variableAttributeLoader.h"

// Constructor
variableAttributeLoader::variableAttributeLoader(){
    setting_primary_field_attributes = true;
    loadVariableAttributes();
    setting_primary_field_attributes = false;
    loadPostProcessorVariableAttributes();
}

// Methods to set the various variable attributes
void variableAttributeLoader::set_variable_name(unsigned int index, std::string name){
    std::pair<unsigned int, std::string> var_pair;
    var_pair.first = index;
    var_pair.second = name;

    if (setting_primary_field_attributes){
        var_name_list.push_back(var_pair);
    }
    else {
        var_name_list_PP.push_back(var_pair);
    }
}
void variableAttributeLoader::set_variable_type(unsigned int index, fieldType var_type){
    std::pair<unsigned int, fieldType> var_pair;
    var_pair.first = index;
    var_pair.second = var_type;

    if (setting_primary_field_attributes){
        var_type_list.push_back(var_pair);
    }
    else {
        var_type_list_PP.push_back(var_pair);
    }
}
void variableAttributeLoader::set_variable_equation_type(unsigned int index, PDEType var_eq_type){
    std::pair<unsigned int, PDEType> var_pair;
    var_pair.first = index;
    var_pair.second = var_eq_type;
    var_eq_type_list.push_back(var_pair);
}
void variableAttributeLoader::set_need_value(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    need_value_list.push_back(var_pair);
}
void variableAttributeLoader::set_need_gradient(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    need_gradient_list.push_back(var_pair);
}
void variableAttributeLoader::set_need_hessian(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    need_hessian_list.push_back(var_pair);
}
void variableAttributeLoader::set_need_value_residual_term(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;

    if (setting_primary_field_attributes){
        need_value_residual_list.push_back(var_pair);
    }
    else {
        need_value_residual_list_PP.push_back(var_pair);
    }
}
void variableAttributeLoader::set_need_gradient_residual_term(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;

    if (setting_primary_field_attributes){
        need_gradient_residual_list.push_back(var_pair);
    }
    else {
        need_gradient_residual_list_PP.push_back(var_pair);
    }
}
void variableAttributeLoader::set_need_value_LHS(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    need_value_list_LHS.push_back(var_pair);
}
void variableAttributeLoader::set_need_gradient_LHS(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    need_gradient_list_LHS.push_back(var_pair);
}
void variableAttributeLoader::set_need_hessian_LHS(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    need_hessian_list_LHS.push_back(var_pair);
}
void variableAttributeLoader::set_need_value_residual_term_LHS(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    need_value_residual_list_LHS.push_back(var_pair);
}
void variableAttributeLoader::set_need_gradient_residual_term_LHS(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    need_gradient_residual_list_LHS.push_back(var_pair);
}
void variableAttributeLoader::set_need_value_PP(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    need_value_list_PP.push_back(var_pair);
}
void variableAttributeLoader::set_need_gradient_PP(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    need_gradient_list_PP.push_back(var_pair);
}
void variableAttributeLoader::set_need_hessian_PP(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    need_hessian_list_PP.push_back(var_pair);
}
void variableAttributeLoader::set_need_value_nucleation(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    need_value_list_nucleation.push_back(var_pair);
}
void variableAttributeLoader::set_allowed_to_nucleate(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    nucleating_variable_list.push_back(var_pair);
}
void variableAttributeLoader::set_output_integral(unsigned int index, bool flag){
    std::pair<unsigned int, bool> var_pair;
    var_pair.first = index;
    var_pair.second = flag;
    output_integral_list.push_back(var_pair);
}
