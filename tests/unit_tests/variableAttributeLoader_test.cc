void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"n");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,PARABOLIC);

	set_need_value					(0,true);
	set_need_gradient				(0,true);
	set_need_hessian					(0,false);

	set_need_value_residual_term		(0,true);
	set_need_gradient_residual_term	(0,true);
}

void variableAttributeLoader::loadPostProcessorVariableAttributes(){
    
}
