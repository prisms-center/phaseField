// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){
	/*
	// Variable 0
	set_variable_name				(0,"f_tot");
	set_variable_type				(0,SCALAR);

	set_need_value_residual_term	(0,true);
	set_need_gradient_residual_term	(0,false);

    set_output_integral         	(0,true);

	// Variable 1
	set_variable_name				(1,"f_el");
	set_variable_type				(1,SCALAR);

	set_need_value_residual_term	(1,true);
	set_need_gradient_residual_term	(1,false);

    set_output_integral         	(1,true);

	// Variable 2
	set_variable_name				(2,"von_mises_stress");
	set_variable_type				(2,SCALAR);

	set_need_value_residual_term	(2,true);
	set_need_gradient_residual_term	(2,false);

	set_output_integral         	(2,true);
	*/
}

// =================================================================================

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
	variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
	const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

		/*
// The order parameter and its derivatives (names here should match those in the macros above)

// The concentration and its derivatives (names here should match those in the macros above)
scalarvalueType c = variable_list.get_scalar_value(0);

// The first order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n1 = variable_list.get_scalar_value(1);
scalargradType n1x = variable_list.get_scalar_gradient(1);

// The second order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n2 = variable_list.get_scalar_value(2);
scalargradType n2x = variable_list.get_scalar_gradient(2);

// The third order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n3 = variable_list.get_scalar_value(3);
scalargradType n3x = variable_list.get_scalar_gradient(3);

// The derivative of the displacement vector (names here should match those in the macros above)
vectorgradType ux = variable_list.get_vector_gradient(4);

scalarvalueType sum_hpV = h1V+h2V+h3V;
scalarvalueType c_alpha = ((B2*c+0.5*(B1-A1)*sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV));
scalarvalueType c_beta  = ((A2*c+0.5*(A1-B1)*(1.0-sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV)));

// Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
// Note: this section can be optimized to reduce recalculations
scalarvalueType cbn1V, cbn2V, cbn3V, cbcV, cbcn1V, cbcn2V, cbcn3V, cacV;

cbcV = faccV/( (constV(1.0)-sum_hpV)*fbccV + (sum_hpV)*faccV );
cacV = fbccV/( (constV(1.0)-sum_hpV)*fbccV + (sum_hpV)*faccV );

cbn1V = hn1V * (c_alpha - c_beta) * cbcV;
cbn2V = hn2V * (c_alpha - c_beta) * cbcV;
cbn3V = hn3V * (c_alpha - c_beta) * cbcV;

cbcn1V = (faccV * (fbccV-faccV) * hn1V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic
cbcn2V = (faccV * (fbccV-faccV) * hn2V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic
cbcn3V = (faccV * (fbccV-faccV) * hn3V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic


scalarvalueType f_chem = (constV(1.0)-(h1V+h2V+h3V))*faV + (h1V+h2V+h3V)*fbV + W*(fbarrierV);

scalarvalueType f_grad = constV(0.0);

for (int i=0; i<dim; i++){
  for (int j=0; j<dim; j++){
	  f_grad += constV(0.5*Kn1[i][j])*n1x[i]*n1x[j] + constV(0.5*Kn2[i][j])*n2x[i]*n2x[j] + constV(0.5*Kn3[i][j])*n3x[i]*n3x[j];
  }
}

// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts1, sfts1c, sfts1cc, sfts2, sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc;

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	// Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c_beta + b_p
	sfts1[i][j] = constV(sfts_linear1[i][j])*c_beta + constV(sfts_const1[i][j]);
	sfts1c[i][j] = constV(sfts_linear1[i][j]) * cbcV;
	sfts1cc[i][j] = constV(0.0);

	sfts2[i][j] = constV(sfts_linear2[i][j])*c_beta + constV(sfts_const2[i][j]);
	sfts2c[i][j] = constV(sfts_linear2[i][j]) * cbcV;
	sfts2cc[i][j] = constV(0.0);

	sfts3[i][j] = constV(sfts_linear3[i][j])*c_beta + constV(sfts_const3[i][j]);
	sfts3c[i][j] = constV(sfts_linear3[i][j]) * cbcV;
	sfts3cc[i][j] = constV(0.0);

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

scalarvalueType total_energy_density = f_chem + f_grad + f_el;


dealii::VectorizedArray<double> vm_stress;
if (dim == 3){
	vm_stress = (S[0][0]-S[1][1])*(S[0][0]-S[1][1]) + (S[1][1]-S[2][2])*(S[1][1]-S[2][2]) + (S[2][2]-S[0][0])*(S[2][2]-S[0][0]);
	vm_stress += constV(6.0)*(S[0][1]*S[0][1] + S[1][2]*S[1][2] + S[2][0]*S[2][0]);
	vm_stress *= constV(0.5);
	vm_stress = std::sqrt(vm_stress);
}
else {
	vm_stress = S[0][0]*S[0][0] - S[0][0]*S[1][1] + S[1][1]*S[1][1] + constV(3.0)*S[0][1]*S[0][1];
	vm_stress = std::sqrt(vm_stress);
}


// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
pp_variable_list.set_scalar_value_residual_term(0, total_energy_density);
pp_variable_list.set_scalar_value_residual_term(1, f_el);
pp_variable_list.set_scalar_value_residual_term(2, vm_stress);
*/
}
