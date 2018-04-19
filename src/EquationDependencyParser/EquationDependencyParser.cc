#include "../../include/EquationDependencyParser.h"


EquationDependencyParser::EquationDependencyParser(
    std::vector<std::string> var_name,
    std::vector<PDEType> var_eq_type,
    std::vector<std::string> sorted_dependencies_value_RHS,
    std::vector<std::string> sorted_dependencies_gradient_RHS,
    std::vector<std::string> sorted_dependencies_value_LHS,
    std::vector<std::string> sorted_dependencies_gradient_LHS)
{
    // Initialize the calculation needed flags to false
    for (unsigned int i=0; i<var_name.size(); i++){
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
    for (unsigned int i=0; i<var_name.size(); i++){

        // First strip excess whitespace
        for(unsigned int j=0; j<sorted_dependencies_value_RHS.at(i).length(); j++){
            if(sorted_dependencies_value_RHS.at(i)[j] == ' ') sorted_dependencies_value_RHS.at(i).erase(j,1);
        }
        // Now check for each variable_eq_type
        if (var_eq_type[i] == PARABOLIC){

            parseDependencyListRHS(var_name, sorted_dependencies_value_RHS.at(i), sorted_dependencies_gradient_RHS.at(i), need_value_explicit_RHS, need_gradient_explicit_RHS, need_hessian_explicit_RHS);

        }
        else if (var_eq_type[i] == ELLIPTIC){

            parseDependencyListRHS(var_name, sorted_dependencies_value_RHS.at(i), sorted_dependencies_gradient_RHS.at(i), need_value_nonexplicit_RHS, need_gradient_nonexplicit_RHS, need_hessian_nonexplicit_RHS);

            parseDependencyListLHS(var_name, sorted_dependencies_value_LHS.at(i), sorted_dependencies_gradient_LHS.at(i), need_value_nonexplicit_LHS, need_gradient_nonexplicit_LHS, need_hessian_nonexplicit_LHS, need_value_change_nonexplicit_LHS, need_gradient_change_nonexplicit_LHS, need_hessian_change_nonexplicit_LHS);

        }
    }
}

void EquationDependencyParser::parseDependencyListRHS(std::vector<std::string> var_name, std::string value_dependencies, std::string gradient_dependencies, std::vector<bool> & need_value, std::vector<bool> & need_gradient, std::vector<bool> & need_hessian){
    std::vector<std::string> split_dependency_list = dealii::Utilities::split_string_list(value_dependencies);
    std::vector<std::string> temp_list = dealii::Utilities::split_string_list(gradient_dependencies);
    split_dependency_list.insert(split_dependency_list.end(),temp_list.begin(),temp_list.end());

    for (unsigned int dep=0; dep<split_dependency_list.size(); dep++){
        bool dependency_entry_assigned = false;

        for (unsigned int var=0; var<var_name.size(); var++){

            std::string grad_var_name = {"grad()"};
            grad_var_name.insert(--grad_var_name.end(),var_name.at(var).begin(),var_name.at(var).end());

            std::string hess_var_name = {"hess()"};
            hess_var_name.insert(--hess_var_name.end(),var_name.at(var).begin(),var_name.at(var).end());

            if (split_dependency_list.at(dep) == var_name.at(var)){
                need_value.at(var) = true;
                dependency_entry_assigned = true;
            }
            else if (split_dependency_list.at(dep) == grad_var_name){
                need_gradient.at(var) = true;
                dependency_entry_assigned = true;
            }
            else if (split_dependency_list.at(dep) == hess_var_name){
                need_hessian.at(var) = true;
                dependency_entry_assigned = true;
            }
        }
        if (!dependency_entry_assigned) {
            std::cerr << "PRISMS-PF Error: Dependency entry " << split_dependency_list.at(dep) << " is not valid." << std::endl;
            abort();
        }
    }
}

void EquationDependencyParser::parseDependencyListLHS(std::vector<std::string> var_name, std::string value_dependencies, std::string gradient_dependencies, std::vector<bool> & need_value, std::vector<bool> & need_gradient, std::vector<bool> & need_hessian, std::vector<bool> & need_value_change, std::vector<bool> & need_gradient_change, std::vector<bool> & need_hessian_change){
    std::vector<std::string> split_dependency_list = dealii::Utilities::split_string_list(value_dependencies);
    std::vector<std::string> temp_list = dealii::Utilities::split_string_list(gradient_dependencies);
    split_dependency_list.insert(split_dependency_list.end(),temp_list.begin(),temp_list.end());

    for (unsigned int dep=0; dep<split_dependency_list.size(); dep++){
        bool dependency_entry_assigned = false;

        for (unsigned int var=0; var<var_name.size(); var++){

            std::string grad_var_name = {"grad()"};
            grad_var_name.insert(--grad_var_name.end(),var_name.at(var).begin(),var_name.at(var).end());

            std::string hess_var_name = {"hess()"};
            hess_var_name.insert(--hess_var_name.end(),var_name.at(var).begin(),var_name.at(var).end());

            std::string val_change_var_name = {"change()"};
            val_change_var_name.insert(--val_change_var_name.end(),var_name.at(var).begin(),var_name.at(var).end());

            std::string grad_change_var_name = {"grad(change())"};
            grad_change_var_name.insert(--(--grad_change_var_name.end()),var_name.at(var).begin(),var_name.at(var).end());

            std::string hess_change_var_name = {"hess(change())"};
            hess_change_var_name.insert(--(--hess_change_var_name.end()),var_name.at(var).begin(),var_name.at(var).end());

            if (split_dependency_list.at(dep) == var_name.at(var)){
                need_value.at(var) = true;
                dependency_entry_assigned = true;
            }
            else if (split_dependency_list.at(dep) == grad_var_name){
                need_gradient.at(var) = true;
                dependency_entry_assigned = true;
            }
            else if (split_dependency_list.at(dep) == hess_var_name){
                need_hessian.at(var) = true;
                dependency_entry_assigned = true;
            }
            else if (split_dependency_list.at(dep) == val_change_var_name){
                need_value_change.at(var) = true;
                dependency_entry_assigned = true;
            }
            else if (split_dependency_list.at(dep) == grad_change_var_name){
                need_gradient_change.at(var) = true;
                dependency_entry_assigned = true;
            }
            else if (split_dependency_list.at(dep) == hess_change_var_name){
                need_hessian_change.at(var) = true;
                dependency_entry_assigned = true;
            }
        }
        if (!dependency_entry_assigned) {
            std::cerr << "PRISMS-PF Error: Dependency entry " << split_dependency_list.at(dep) << " is not valid." << std::endl;
            abort();
        }
    }
}
