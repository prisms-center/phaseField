// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){

	// Variable 0
	set_variable_name				(0,"f_tot");
	set_variable_type				(0,SCALAR);

	set_need_value_residual_term	(0,true);
	set_need_gradient_residual_term	(0,false);

    set_output_integral         	(0,true);

}

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

scalarvalueType f_tot = constV(0.0);

//c
scalarvalueType c = variable_list.get_scalar_value(0);

//n1
scalarvalueType n = variable_list.get_scalar_value(1);
scalargradType nx = variable_list.get_scalar_gradient(1);

//biharm
scalarvalueType biharm = variable_list.get_scalar_value(2);

scalarvalueType faV = 0.5*c*c/16.0;
scalarvalueType fbV = 0.5*(c-1.0)*(c-1.0)/16.0;
scalarvalueType hV = 3.0*n*n-2.0*n*n*n;

scalarvalueType normgradn = std::sqrt(nx.norm_square());
scalargradType normal = nx/(normgradn+constV(1.0e-16));

scalarvalueType gamma;
if (dim == 2){
    gamma = 1.0+epsilonM*(4.0*(normal[0]*normal[0]*normal[0]*normal[0]+normal[1]*normal[1]*normal[1]*normal[1])-3.0);
}
else {
    gamma = 1.0+epsilonM*(4.0*(normal[0]*normal[0]*normal[0]*normal[0]+normal[1]*normal[1]*normal[1]*normal[1]+normal[2]*normal[2]*normal[2]*normal[2])-3.0);
}

scalarvalueType f_chem = (constV(1.0)-hV)*faV + hV*fbV;

// anisotropy code
scalarvalueType f_grad = constV(0.5)*gamma*gamma*nx*nx;

scalarvalueType f_reg = constV(0.5*delta2)*biharm*biharm;

f_tot = f_chem + f_grad + f_reg;

// end anisotropy code
f_tot = f_chem + f_grad;

// Residuals for the equation to evolve the order parameter 
pp_variable_list.set_scalar_value_residual_term(0, f_tot);


}
