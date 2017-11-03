// List of variables and residual equations for the Precipitate Evolution example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"c");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,PARABOLIC);

	set_need_value					(0,true);
	set_need_gradient				(0,true);
	set_need_hessian				(0,false);

	set_need_value_residual_term	(0,true);
	set_need_gradient_residual_term	(0,true);

	set_need_value_LHS				(0,false);
	set_need_gradient_LHS			(0,false);
	set_need_hessian_LHS			(0,false);
	set_need_value_residual_term_LHS	(0,false);
	set_need_gradient_residual_term_LHS	(0,false);

	set_allowed_to_nucleate			(0, false);
	set_need_value_nucleation		(0, true);

	// Variable 1
	set_variable_name				(1,"n1");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,PARABOLIC);

	set_need_value					(1,true);
	set_need_gradient				(1,true);
	set_need_hessian				(1,false);

	set_need_value_residual_term	(1,true);
	set_need_gradient_residual_term	(1,true);

	set_need_value_LHS				(1,true);
	set_need_gradient_LHS			(1,false);
	set_need_hessian_LHS			(1,false);
	set_need_value_residual_term_LHS	(1,false);
	set_need_gradient_residual_term_LHS	(1,false);

	set_allowed_to_nucleate			(1, true);
	set_need_value_nucleation		(1, true);

	// Variable 2
	set_variable_name				(2,"n2");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,PARABOLIC);

	set_need_value					(2,true);
	set_need_gradient				(2,true);
	set_need_hessian				(2,false);

	set_need_value_residual_term	(2,true);
	set_need_gradient_residual_term	(2,true);

	set_need_value_LHS				(2,true);
	set_need_gradient_LHS			(2,false);
	set_need_hessian_LHS			(2,false);
	set_need_value_residual_term_LHS	(2,false);
	set_need_gradient_residual_term_LHS	(2,false);

	set_allowed_to_nucleate			(2, true);
	set_need_value_nucleation		(2, true);

	// Variable 3
	set_variable_name				(3,"n3");
	set_variable_type				(3,SCALAR);
	set_variable_equation_type		(3,PARABOLIC);

	set_need_value					(3,true);
	set_need_gradient				(3,true);
	set_need_hessian				(3,false);

	set_need_value_residual_term	(3,true);
	set_need_gradient_residual_term	(3,true);

	set_need_value_LHS				(3,true);
	set_need_gradient_LHS			(3,false);
	set_need_hessian_LHS			(3,false);
	set_need_value_residual_term_LHS	(3,false);
	set_need_gradient_residual_term_LHS	(3,false);

	set_allowed_to_nucleate			(3, true);
	set_need_value_nucleation		(3, true);

	// Variable 4
	set_variable_name				(4,"u");
	set_variable_type				(4,VECTOR);
	set_variable_equation_type		(4,ELLIPTIC);

	set_need_value					(4,false);
	set_need_gradient				(4,true);
	set_need_hessian				(4,false);

	set_need_value_residual_term	(4,false);
	set_need_gradient_residual_term	(4,true);

	set_need_value_LHS				(4,false);
	set_need_gradient_LHS			(4,true);
	set_need_hessian_LHS			(4,false);
	set_need_value_residual_term_LHS	(4,false);
	set_need_gradient_residual_term_LHS	(4,true);

}

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".

#define faV (A2*c_alpha*c_alpha + A1*c_alpha + A0)
#define facV (2.0*A2*c_alpha +A1)
#define faccV (2.0*A2)
#define fbV (B2*c_beta*c_beta + B1*c_beta + B0)
#define fbcV (2.0*B2*c_beta +B1)
#define fbccV (2.0*B2)

#define h1V (3.0*n1*n1-2.0*n1*n1*n1)
#define hn1V (6.0*n1-6.0*n1*n1)
#define h2V (3.0*n2*n2-2.0*n2*n2*n2)
#define hn2V (6.0*n2-6.0*n2*n2)
#define h3V (3.0*n3*n3-2.0*n3*n3*n3)
#define hn3V (6.0*n3-6.0*n3*n3)

// This double-well function can be used to tune the interfacial energy
#define fbarrierV (n1*n1-2.0*n1*n1*n1+n1*n1*n1*n1 + n2*n2-2.0*n2*n2*n2+n2*n2*n2*n2 + n3*n3-2.0*n3*n3*n3+n3*n3*n3*n3 + 5.0*(n1*n1*n2*n2 + n1*n1*n3*n3 + n2*n2*n3*n3) + 5.0*n1*n1*n2*n2*n3*n3)
#define fbarriern1V (2.0*n1-6.0*n1*n1+4.0*n1*n1*n1 + 10.0*n1*(n2*n2+n3*n3) + 10.0*n1*n2*n2*n3*n3)
#define fbarriern2V (2.0*n2-6.0*n2*n2+4.0*n2*n2*n2 + 10.0*n2*(n1*n1+n3*n3) + 10.0*n2*n1*n1*n3*n3)
#define fbarriern3V (2.0*n3-6.0*n3*n3+4.0*n3*n3*n3 + 10.0*n3*(n2*n2+n1*n1) + 10.0*n3*n2*n2*n1*n1)

// Residuals
#define rcV   (c)
#define rcxTemp ( cx + (c_alpha-c_beta)*(hn1V*n1x + hn2V*n2x + hn3V*n3x) + grad_mu_el * ((sum_hpV)*faccV+(constV(1.0)-sum_hpV)*fbccV)/constV(faccV*fbccV))
#define rcxV  (constV(-userInputs.dtValue)*McV*rcxTemp)

#define rn1V   (n1-constV(userInputs.dtValue*Mn1V)*gamma*( (fbV-faV)*hn1V - (c_beta-c_alpha)*facV*hn1V + W*fbarriern1V + nDependentMisfitAC1 + heterMechAC1))
#define rn1xV  (constV(-userInputs.dtValue*Mn1V)*gamma*Knx1)

#define rn2V   (n2-constV(userInputs.dtValue*Mn2V)*gamma*( (fbV-faV)*hn2V - (c_beta-c_alpha)*facV*hn2V + W*fbarriern2V + nDependentMisfitAC2 + heterMechAC2))
#define rn2xV  (constV(-userInputs.dtValue*Mn2V)*gamma*Knx2)

#define rn3V   (n3-constV(userInputs.dtValue*Mn3V)*gamma*( (fbV-faV)*hn3V - (c_beta-c_alpha)*facV*hn3V + W*fbarriern3V + nDependentMisfitAC3 + heterMechAC3))
#define rn3xV  (constV(-userInputs.dtValue*Mn3V)*gamma*Knx3)

// =================================================================================
// residualRHS
// =================================================================================
// This function calculates the residual equations for each variable. It takes
// "modelVariablesList" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of
// that quadrature point is given by "q_point_loc". The function outputs
// "modelResidualsList", a list of the value and gradient terms of the residual for
// each residual equation. The index for each variable in these lists corresponds to
// the order it is defined at the top of this file (starting at 0).
template <int dim, int degree>
void customPDE<dim,degree>::residualRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {


// The concentration and its derivatives (names here should match those in the macros above)
scalarvalueType c = variable_list.get_scalar_value(0);
scalargradType cx = variable_list.get_scalar_gradient(0);

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
vectorgradType ruxV;

vectorhessType uxx;

if (c_dependent_misfit == true){
	uxx = variable_list.get_vector_hessian(4);
}

// -------------------------------------------------
// Nucleation expressions
// -------------------------------------------------
std::vector<dealii::VectorizedArray<double> > source_terms(3,constV(0.0));
dealii::VectorizedArray<double> gamma = constV(1.0);
seedNucleus(q_point_loc,source_terms,gamma);

// -------------------------------------------------

// Calculate c_alpha and c_beta from c
scalarvalueType sum_hpV = h1V+h2V+h3V;
scalarvalueType c_alpha = ((B2*c+0.5*(B1-A1)*sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV));
scalarvalueType c_beta  = ((A2*c+0.5*(A1-B1)*(1.0-sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV)));

// Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
// Note: this section can be optimized to reduce recalculations (note: that doesn't change the performance)
scalarvalueType cbn1V, cbn2V, cbn3V, cbcV, cbcn1V, cbcn2V, cbcn3V;


cbcV = faccV/( (constV(1.0)-sum_hpV)*fbccV + (sum_hpV)*faccV );

cbn1V = hn1V * (c_alpha - c_beta) * cbcV;
cbn2V = hn2V * (c_alpha - c_beta) * cbcV;
cbn3V = hn3V * (c_alpha - c_beta) * cbcV;

cbcn1V = (faccV * (fbccV-faccV) * hn1V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic
cbcn2V = (faccV * (fbccV-faccV) * hn2V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic
cbcn3V = (faccV * (fbccV-faccV) * hn3V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic


// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts1, sfts1c, sfts1cc, sfts2, sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc;

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	// Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c_beta + b_p
	sfts1[i][j] = constV(sfts_linear1[i][j])*c_beta + constV(sfts_const1[i][j]);
	//sfts1c[i][j] = constV(sfts_linear1[i][j]) * cbcV;
	//sfts1cc[i][j] = constV(0.0);

	sfts2[i][j] = constV(sfts_linear2[i][j])*c_beta + constV(sfts_const2[i][j]);
	//sfts2c[i][j] = constV(sfts_linear2[i][j]) * cbcV;
	//sfts2cc[i][j] = constV(0.0);

	sfts3[i][j] = constV(sfts_linear3[i][j])*c_beta + constV(sfts_const3[i][j]);
	//sfts3c[i][j] = constV(sfts_linear3[i][j]) * cbcV;
	//sfts3cc[i][j] = constV(0.0);

}
}

//compute E2=(E-E0)
dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V + sfts2[i][j]*h2V + sfts3[i][j]*h3V );

}
}

//compute stress
//S=C*(E-E0)
// Compute stress tensor (which is equal to the residual, Rux)
dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

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


// Fill residual corresponding to mechanics
// R=-C*(E-E0)

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  ruxV[i][j] = - S[i][j];
}
}

// Compute one of the stress terms in the order parameter chemical potential, nDependentMisfitACp = -C*(E-E0)*(E0_n)
dealii::VectorizedArray<double> nDependentMisfitAC1=constV(0.0);
dealii::VectorizedArray<double> nDependentMisfitAC2=constV(0.0);
dealii::VectorizedArray<double> nDependentMisfitAC3=constV(0.0);

dealii::VectorizedArray<double> E4[dim][dim]; // Intermediate variable

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	E4[i][j] = constV(sfts_linear1[i][j])*h1V + constV(sfts_linear2[i][j])*h2V + constV(sfts_linear3[i][j])*h3V;

	nDependentMisfitAC1 -= S[i][j]*(cbn1V*E4[i][j] + sfts1[i][j]*hn1V);
	nDependentMisfitAC2 -= S[i][j]*(cbn2V*E4[i][j] + sfts2[i][j]*hn2V);
	nDependentMisfitAC3 -= S[i][j]*(cbn3V*E4[i][j] + sfts3[i][j]*hn3V);
}
}

// Compute the other stress term in the order parameter chemical potential, heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
dealii::VectorizedArray<double> heterMechAC1=constV(0.0);
dealii::VectorizedArray<double> heterMechAC2=constV(0.0);
dealii::VectorizedArray<double> heterMechAC3=constV(0.0);
dealii::VectorizedArray<double> S2[dim][dim];

if (n_dependent_stiffness == true){
	computeStress<dim>(CIJ_Beta1-CIJ_Mg, E2, S2);
	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			heterMechAC1 += S2[i][j]*E2[i][j];
		}
	}
	heterMechAC1 *= 0.5*hn1V;

	computeStress<dim>(CIJ_Beta2-CIJ_Mg, E2, S2);
	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			heterMechAC2 += S2[i][j]*E2[i][j];
		}
	}
	heterMechAC2 *= 0.5*hn2V;

	computeStress<dim>(CIJ_Beta3-CIJ_Mg, E2, S2);
	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			heterMechAC3 += S2[i][j]*E2[i][j];
		}
	}
	heterMechAC3 *= 0.5*hn3V;
}

// compute the stress term in the gradient of the concentration chemical potential, grad_mu_el = [C*(E-E0)*E0c]x, must be a vector with length dim
scalargradType grad_mu_el;

if (c_dependent_misfit == true){
	dealii::VectorizedArray<double> E3[dim][dim], S3[dim][dim]; // Intermediate variables

	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			E3[i][j] =  -( sfts1c[i][j]*h1V + sfts2c[i][j]*h2V + sfts3c[i][j]*h3V );
		}
	}

	if (n_dependent_stiffness == true){
		computeStress<dim>(CIJ_combined, E3, S3);
	}
	else{
		computeStress<dim>(CIJ_Mg, E3, S3);
	}

	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			for (unsigned int k=0; k<dim; k++){
				grad_mu_el[k] += S3[i][j] * (constV(0.5)*(uxx[i][j][k]+uxx[j][i][k]) + E3[i][j]*cx[k]
										- ( (sfts1[i][j]*hn1V + cbn1V*E4[i][j])*n1x[k] + (sfts2[i][j]*hn2V + cbn2V*E4[i][j])*n2x[k]
										+ (sfts3[i][j]*hn3V + cbn3V*E4[i][j])*n3x[k]) );

				grad_mu_el[k]+= - S[i][j] * ( (sfts1c[i][j]*hn1V + cbcn1V*E4[i][j])*n1x[k]
										+ (sfts2c[i][j]*hn2V+ cbcn2V*E4[i][j])*n2x[k]
										+ (sfts3c[i][j]*hn3V+ cbcn3V*E4[i][j])*n3x[k]
										+ ( sfts1cc[i][j]*h1V + sfts2cc[i][j]*h2V + sfts3cc[i][j]*h3V )*cx[k]);

				if (n_dependent_stiffness == true){
					computeStress<dim>(CIJ_Beta1-CIJ_Mg, E2, S2);
					grad_mu_el[k]+= - S2[i][j] * (cbcV*E4[i][j]*hn1V*n1x[k]);
					computeStress<dim>(CIJ_Beta2-CIJ_Mg, E2, S2);
					grad_mu_el[k]+= - S2[i][j] * (cbcV*E4[i][j]*hn2V*n2x[k]);
					computeStress<dim>(CIJ_Beta3-CIJ_Mg, E2, S2);
					grad_mu_el[k]+= - S2[i][j] * (cbcV*E4[i][j]*hn3V*n3x[k]);

				}
			}
		}
	}
}


//compute K*nx
scalargradType Knx1, Knx2, Knx3;
for (unsigned int a=0; a<dim; a++) {
Knx1[a]=0.0;
Knx2[a]=0.0;
Knx3[a]=0.0;
for (unsigned int b=0; b<dim; b++){
	  Knx1[a]+=constV(Kn1[a][b])*n1x[b];
	  Knx2[a]+=constV(Kn2[a][b])*n2x[b];
	  Knx3[a]+=constV(Kn3[a][b])*n3x[b];
}
}

variable_list.set_scalar_value_residual_term(0,rcV);
variable_list.set_scalar_gradient_residual_term(0,rcxV);

variable_list.set_scalar_value_residual_term(1,rn1V + source_terms[0]);
variable_list.set_scalar_gradient_residual_term(1,rn1xV);

variable_list.set_scalar_value_residual_term(2,rn2V + source_terms[1]);
variable_list.set_scalar_gradient_residual_term(2,rn2xV);

variable_list.set_scalar_value_residual_term(3,rn3V + source_terms[2]);
variable_list.set_scalar_gradient_residual_term(3,rn3xV);

variable_list.set_vector_gradient_residual_term(4,ruxV);

}

// =================================================================================
// seedNucleus: a function particular to this app
// =================================================================================
template <int dim,int degree>
void customPDE<dim,degree>::seedNucleus(const dealii::Point<dim, dealii::VectorizedArray<double> > & q_point_loc,
	std::vector<dealii::VectorizedArray<double> > & source_terms,
	dealii::VectorizedArray<double> & gamma) const {

		for (typename std::vector<nucleus<dim> >::const_iterator thisNucleus=this->nuclei.begin(); thisNucleus!=this->nuclei.end(); ++thisNucleus){

			if (thisNucleus->seededTime + thisNucleus->seedingTime > this->currentTime){

				// Get the rotation angle
				double theta;
				if (thisNucleus->orderParameterIndex == 1){
					theta = -pi/3.0;
				}
				else if (thisNucleus->orderParameterIndex == 2){
					theta = 0.0;
				}
				else {
					theta = pi/3.0;
				}

				// Calculate the weighted distance function to the order parameter freeze boundary (weighted_dist = 1.0 on that boundary)
				// Allow rotation of the ellipsoid in the x-y plane depending on the order parameter (as given by theta)
				dealii::VectorizedArray<double> weighted_dist = constV(0.0);
				for (unsigned int i=0; i<dim; i++){
					dealii::VectorizedArray<double> temp = (thisNucleus->center(i) - q_point_loc(i));

					if (userInputs.BC_list[1].var_BC_type[2*i]==PERIODIC){
						for (unsigned j=0; j<gamma.n_array_elements;j++)
						temp[j] -= round(temp[j]/userInputs.domain_size[i])*userInputs.domain_size[i];
					}

					if (i == 0){
						dealii::VectorizedArray<double> temp_y = (thisNucleus->center(1) - q_point_loc(1));
						temp = temp*std::cos(theta)+temp_y*std::sin(theta);
					}
					else if (i == 1){
						dealii::VectorizedArray<double> temp_x = (thisNucleus->center(0) - q_point_loc(0));
						temp = -temp_x*std::sin(theta)+temp*std::cos(theta);
					}

					temp=temp/userInputs.order_parameter_freeze_semiaxes[i];
					weighted_dist += temp*temp;
				}

				for (unsigned i=0; i<gamma.n_array_elements;i++){
					if (weighted_dist[i] <= 1.0){
						gamma[i] = 0.0;

						// Seed a nucleus if it was added to the list of nuclei this time step
						if (thisNucleus->seedingTimestep == this->currentIncrement){

							// Find the weighted distance to the outer edge of the nucleus and use it to calculate the order parameter source term
							// Allow rotation of the ellipsoid in the x-y plane depending on the order parameter (as given by theta)
							double r = 0.0;
							double avg_semiaxis = 0.0;
							for (unsigned int j=0; j<dim; j++){
								double temp = (thisNucleus->center(j) - q_point_loc(j)[i]);
								if (userInputs.BC_list[1].var_BC_type[2*i]==PERIODIC){
									double domsize_j = userInputs.domain_size[j];
									temp=temp-round(temp/domsize_j)*domsize_j;
								}

								if (j == 0){
									double temp_y = (thisNucleus->center(1) - q_point_loc(1)[i]);
									temp = temp*std::cos(theta)+temp_y*std::sin(theta);
								}
								else if (j == 1){
									double temp_x = (thisNucleus->center(0) - q_point_loc(0)[i]);
									temp = -temp_x*std::sin(theta)+temp*std::cos(theta);
								}

								temp=temp/thisNucleus->semiaxes[j];
								r += temp*temp;
								avg_semiaxis += thisNucleus->semiaxes[j];
							}

							r = sqrt(r);
							avg_semiaxis /= dim;
							if (thisNucleus->orderParameterIndex == 1){
								source_terms[0][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/interface_coeff));
							}
							else if (thisNucleus->orderParameterIndex == 2){
								source_terms[1][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/interface_coeff));
							}
							else {
								source_terms[2][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/interface_coeff));
							}
						}
					}
				}
			}
		}
	}

// =================================================================================
// residualLHS (needed only if at least one equation is elliptic)
// =================================================================================
// This function calculates the residual equations for the iterative solver for
// elliptic equations.for each variable. It takes "modelVariablesList" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is given
// by "q_point_loc". The function outputs "modelRes", the value and gradient terms of
// for the left-hand-side of the residual equation for the iterative solver. The
// index for each variable in these lists corresponds to the order it is defined at
// the top of this file (starting at 0), not counting variables that have
// "need_val_LHS", "need_grad_LHS", and "need_hess_LHS" all set to "false". If there
// are multiple elliptic equations, conditional statements should be used to ensure
// that the correct residual is being submitted. The index of the field being solved
// can be accessed by "this->currentFieldIndex".
template <int dim, int degree>
void customPDE<dim,degree>::residualLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {


//n1
scalarvalueType n1 = variable_list.get_scalar_value(1);

//n2
scalarvalueType n2 = variable_list.get_scalar_value(2);

//n3
scalarvalueType n3 = variable_list.get_scalar_value(3);

vectorgradType ruxV;

// Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the dealii "symmetrize" function
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > E = symmetrize(variable_list.get_vector_gradient(4));

// Compute stress tensor (which is equal to the residual, Rux)
if (n_dependent_stiffness == true){
	dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > CIJ_combined;
	CIJ_combined = CIJ_Mg*(constV(1.0)-h1V-h2V-h3V) + CIJ_Beta1*(h1V) + CIJ_Beta2*(h2V) + CIJ_Beta3*(h3V);

	computeStress<dim>(CIJ_combined, E, ruxV);
}
else{
	computeStress<dim>(CIJ_Mg, E,ruxV);
}

variable_list.set_vector_gradient_residual_term(4,ruxV);

}
