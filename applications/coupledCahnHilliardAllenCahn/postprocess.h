// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){

	// Variable 0
	set_variable_name				(0,"f_tot");
	set_variable_type				(0,SCALAR);

    set_dependencies_value_residual_term_RHS(0, "c,n,grad(n)");
    set_dependencies_gradient_residual_term_RHS(0, "");

    set_output_integral         	(0,true);

}

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//c
scalarvalueType c = variable_list.get_scalar_value(0);

//n
scalarvalueType n = variable_list.get_scalar_value(1);
scalargradType nx = variable_list.get_scalar_gradient(1);

// Free energy for each phase
scalarvalueType fa = constV(2.0)*c*c;
scalarvalueType fb = constV(2.0)*(c*c - 2.0*c + constV(1.0));

// Interpolation function
scalarvalueType h = (3.0*n*n-2.0*n*n*n);

// The homogenous free energy
scalarvalueType f_chem = (constV(1.0)-h)*fa + h*fb;

// The gradient free energy
scalarvalueType f_grad = constV(0.5*Kn)*nx*nx;

// The total free energy
scalarvalueType f_tot;
f_tot = f_chem + f_grad;

// Residuals for the equation to evolve the order parameter
pp_variable_list.set_scalar_value_residual_term(0, f_tot);


}
