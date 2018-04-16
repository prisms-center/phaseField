// Class to hold the variable attributes that will be passed to a userInputParameters object
#ifndef VARIABLEATTRIBUTELOADER_H
#define VARIABLEATTRIBUTELOADER_H

#include "dealIIheaders.h"
#include "model_variables.h"
#include "varTypeEnums.h"

class variableAttributeLoader
{
public:
    // Constructor
    variableAttributeLoader();

    // Methods where the attributes are set
    void loadVariableAttributes();
    void loadPostProcessorVariableAttributes();

    // Methods to set the parameter_attributes
	bool setting_primary_field_attributes;
	void set_variable_name(unsigned int index, std::string name);
	void set_variable_type(unsigned int index, fieldType);
	void set_variable_equation_type(unsigned int index, PDEType);

	void set_need_value(unsigned int index, bool);
	void set_need_gradient(unsigned int index, bool);
	void set_need_hessian(unsigned int index, bool);
	void set_need_value_residual_term(unsigned int index, bool);
	void set_need_gradient_residual_term(unsigned int index, bool);

	void set_need_value_LHS(unsigned int index, bool);
	void set_need_gradient_LHS(unsigned int index, bool);
	void set_need_hessian_LHS(unsigned int index, bool);
	void set_need_value_residual_term_LHS(unsigned int index, bool);
	void set_need_gradient_residual_term_LHS(unsigned int index, bool);

    void set_need_value_change_LHS(unsigned int index, bool);
	void set_need_gradient_change_LHS(unsigned int index, bool);
	void set_need_hessian_change_LHS(unsigned int index, bool);

    void set_need_value_PP(unsigned int index, bool);
	void set_need_gradient_PP(unsigned int index, bool);
	void set_need_hessian_PP(unsigned int index, bool);

    void set_need_value_nucleation(unsigned int index, bool);
    void set_allowed_to_nucleate(unsigned int index, bool);

    void set_output_integral(unsigned int index, bool);

    void set_equations_are_nonlinear(bool);

    // Variable inputs (v2.0)
    std::vector<std::pair<unsigned int, std::string> > var_name_list;
    std::vector<std::pair<unsigned int, fieldType> > var_type_list;
    std::vector<std::pair<unsigned int, PDEType> > var_eq_type_list;
    std::vector<std::pair<unsigned int, bool> > need_value_list;
    std::vector<std::pair<unsigned int, bool> > need_gradient_list;
    std::vector<std::pair<unsigned int, bool> > need_hessian_list;
    std::vector<std::pair<unsigned int, bool> > need_value_residual_list;
    std::vector<std::pair<unsigned int, bool> > need_gradient_residual_list;
    std::vector<std::pair<unsigned int, bool> > need_value_list_LHS;
    std::vector<std::pair<unsigned int, bool> > need_gradient_list_LHS;
    std::vector<std::pair<unsigned int, bool> > need_hessian_list_LHS;
    std::vector<std::pair<unsigned int, bool> > need_value_residual_list_LHS;
    std::vector<std::pair<unsigned int, bool> > need_gradient_residual_list_LHS;

    std::vector<std::pair<unsigned int, bool> > need_value_change_list_LHS;
    std::vector<std::pair<unsigned int, bool> > need_gradient_change_list_LHS;
    std::vector<std::pair<unsigned int, bool> > need_hessian_change_list_LHS;

    std::vector<std::pair<unsigned int, bool> > need_value_list_PP;
    std::vector<std::pair<unsigned int, bool> > need_gradient_list_PP;
    std::vector<std::pair<unsigned int, bool> > need_hessian_list_PP;
    std::vector<std::pair<unsigned int, bool> > need_value_list_nucleation;
    std::vector<std::pair<unsigned int, bool> > nucleating_variable_list;

    std::vector<std::pair<unsigned int, std::string> > var_name_list_PP;
    std::vector<std::pair<unsigned int, fieldType> > var_type_list_PP;
    std::vector<std::pair<unsigned int, bool> > output_integral_list;
    std::vector<std::pair<unsigned int, bool> > need_value_residual_list_PP;
    std::vector<std::pair<unsigned int, bool> > need_gradient_residual_list_PP;

    bool equations_are_nonlinear;

    void set_dependencies_value_RHS(unsigned int index, std::string dependencies);
    void set_dependencies_gradient_RHS(unsigned int index, std::string dependencies);
    void set_dependencies_value_LHS(unsigned int index, std::string dependencies);
    void set_dependencies_gradient_LHS(unsigned int index, std::string dependencies);

    std::vector<std::pair<unsigned int, std::string> > var_eq_dependencies_value_RHS;
    std::vector<std::pair<unsigned int, std::string> > var_eq_dependencies_gradient_RHS;
    std::vector<std::pair<unsigned int, std::string> > var_eq_dependencies_value_LHS;
    std::vector<std::pair<unsigned int, std::string> > var_eq_dependencies_gradient_LHS;

};


#endif
