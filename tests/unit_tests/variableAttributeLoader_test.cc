void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"n");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,PARABOLIC);


    set_dependencies_value_residual_term_RHS(0, "n");
    set_dependencies_gradient_residual_term_RHS(0, "grad(n)");

}

void variableAttributeLoader::loadPostProcessorVariableAttributes(){

}
