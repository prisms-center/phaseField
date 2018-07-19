// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){

	// Variable 1
	set_variable_name				(0,"f_tot");
	set_variable_type				(0,SCALAR);

    set_dependencies_value_residual_term_RHS(0, "u, phi, grad(phi)");
    set_dependencies_gradient_residual_term_RHS(0, "");

    set_output_integral         	(0,true);

}

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// The order parameter and its derivatives
// The temperature and its derivatives
scalarvalueType u = variable_list.get_scalar_value(0);

// The order parameter and its derivatives
scalarvalueType phi = variable_list.get_scalar_value(1);
scalargradType phix = variable_list.get_scalar_gradient(1);

double lambda = (D/0.6267/W0/W0);

scalarvalueType f_tot = constV(0.0);

// The homogenous free energy
scalarvalueType f_chem = -0.5 * phi*phi + 0.25*phi*phi*phi*phi + lambda*u*phi*(1.0-2.0/3.0*phi*phi+1.0/5.0/phi*phi*phi*phi);

// The azimuthal angle
scalarvalueType theta;
for (unsigned i=0; i< phi.n_array_elements;i++){
	theta[i] = std::atan2(phix[1][i],phix[0][i]);
}

scalarvalueType W = constV(W0)*(constV(1.0)+constV(epsilonM)*std::cos(constV(mult)*(theta-constV(theta0))));

// The gradient free energy
scalarvalueType f_grad = constV(0.5)*W*W * (phix[0]*phix[0] + phix[1]*phix[1]);

// The total free energy
f_tot = f_chem + f_grad;


// Residuals for the equation to evolve the order parameter
pp_variable_list.set_scalar_value_residual_term(0, f_tot);


}
