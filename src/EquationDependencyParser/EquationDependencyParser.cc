#include "../../include/EquationDependencyParser.h"


void EquationDependencyParser::parseEquationDependencies(
    std::vector<std::string> sorted_dependencies_value_RHS,
    std::vector<std::string> sorted_dependencies_gradient_RHS,
    std::vector<std::string> sorted_dependencies_value_LHS,
    std::vector<std::string> sorted_dependencies_gradient_LHS,
    std::vector<bool> need_value_explicit_RHS,
    std::vector<bool> need_gradient_explicit_RHS,
    std::vector<bool> need_hessian_explicit_RHS,
    std::vector<bool> need_value_nonexplicit_RHS,
    std::vector<bool> need_gradient_nonexplicit_RHS,
    std::vector<bool> need_hessian_nonexplicit_RHS,
    std::vector<bool> need_value_nonexplicit_LHS,
    std::vector<bool> need_gradient_nonexplicit_LHS,
    std::vector<bool> need_hessian_nonexplicit_LHS,
    std::vector<bool> need_value_change_nonexplicit_LHS,
    std::vector<bool> need_gradient_change_nonexplicit_LHS,
    std::vector<bool> need_hessian_change_nonexplicit_LHS)
{
    for (unsigned int i=0; i<number_of_variables; i++){
        // First strip excess whitespace
        for(unsigned int j=0; j<sorted_dependencies_value_RHS.at(i).length(); j++){
            if(sorted_dependencies_value_RHS.at(i)[j] == ' ') sorted_dependencies_value_RHS.at(i).erase(j,1);
        }
        // Now check for each variable_eq_type


    }



}


SolverToleranceType SolverParametersBase::getToleranceType(unsigned int index){
    return tolerance_type_list.at(getEquationIndex(index));
}
