// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){

	// Variable 0
	set_variable_name				(0,"f_el");
	set_variable_type				(0,SCALAR);

    set_dependencies_value_term_RHS(0, "c, n1, n2, n3, grad(u)");
    set_dependencies_gradient_term_RHS(0, "");

    set_output_integral         	(0,false);

	// Variable 1
	set_variable_name				(1,"e_11");
	set_variable_type				(1,SCALAR);

    set_dependencies_value_term_RHS(1, "c, n1, n2, n3, grad(u)");
    set_dependencies_gradient_term_RHS(1, "");

	set_output_integral         	(1,false);

    // Variable 2
	set_variable_name				(2,"e_21");
	set_variable_type				(2,SCALAR);

    set_dependencies_value_term_RHS(2, "c, n1, n2, n3, grad(u)");
    set_dependencies_gradient_term_RHS(2, "");

	set_output_integral         	(2,false);

    // Variable 3
    set_variable_name				(3,"e_22");
    set_variable_type				(3,SCALAR);

    set_dependencies_value_term_RHS(3, "c, n1, n2, n3, grad(u)");
    set_dependencies_gradient_term_RHS(3, "");

    set_output_integral         	(3,false);





}

// =================================================================================

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
	variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
	const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {


// The order parameter and its derivatives (names here should match those in the macros above)

// The concentration and its derivatives (names here should match those in the macros above)
scalarvalueType c = variable_list.get_scalar_value(0);

// The first order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n1 = variable_list.get_scalar_value(1);

// The second order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n2 = variable_list.get_scalar_value(2);

// The third order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n3 = variable_list.get_scalar_value(3);

// The derivative of the displacement vector (names here should match those in the macros above)
vectorgradType ux = variable_list.get_vector_gradient(4);

scalarvalueType sum_hpV = h1V+h2V+h3V;
scalarvalueType c_alpha = ((B2*c+0.5*(B1-A1)*sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV));
scalarvalueType c_beta  = ((A2*c+0.5*(A1-B1)*(1.0-sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV)));


// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts1, sfts1c, sfts1cc, sfts2, sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc;

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	// Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c_beta + b_p
	sfts1[i][j] = constV(sfts_linear1[i][j])*c_beta + constV(sfts_const1[i][j]);

	sfts2[i][j] = constV(sfts_linear2[i][j])*c_beta + constV(sfts_const2[i][j]);

	sfts3[i][j] = constV(sfts_linear3[i][j])*c_beta + constV(sfts_const3[i][j]);

}
}


//compute E2=(E-E0)
scalarvalueType E2[dim][dim], S[dim][dim];

for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V + sfts2[i][j]*h2V + sfts3[i][j]*h3V);

  }
}

//compute stress
//S=C*(E-E0)
scalarvalueType CIJ_combined[2*dim-1+dim/3][2*dim-1+dim/3];

if (n_dependent_stiffness == true){
for (unsigned int i=0; i<2*dim-1+dim/3; i++){
	  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
		  CIJ_combined[i][j] = CIJ_Mg[i][j]*(constV(1.0)-sum_hpV) + CIJ_Beta1[i][j]*(h1V) + CIJ_Beta2[i][j]*(h2V) + CIJ_Beta3[i][j]*(h3V);
	  }
}
computeStress<dim>(CIJ_combined, E2, S);
}
else{
computeStress<dim>(CIJ_Mg, E2, S);
}

dealii::VectorizedArray<double> f_el = constV(0.0);

for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  f_el += constV(0.5) * S[i][j]*E2[i][j];
  }
}

// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)

pp_variable_list.set_scalar_value_term_RHS(0, f_el);
pp_variable_list.set_scalar_value_term_RHS(1, S[0][0]);
pp_variable_list.set_scalar_value_term_RHS(2, S[1][0]);
pp_variable_list.set_scalar_value_term_RHS(3, S[1][1]);


}
