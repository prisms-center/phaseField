// List of variables and residual equations for the Precipitate Evolution example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"c");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "c");
    set_dependencies_gradient_term_RHS(0, "n1, n2, n3, grad(c), grad(n1), grad(n2), grad(n3), grad(u)");


	// Variable 1
	set_variable_name				(1,"n1");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(1, "c, n1, n2, n3, grad(u)");
    set_dependencies_gradient_term_RHS(1, "grad(n1)");


	// Variable 2
	set_variable_name				(2,"n2");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(2, "c, n1, n2, n3, grad(u)");
    set_dependencies_gradient_term_RHS(2, "grad(n2)");

	// Variable 3
	set_variable_name				(3,"n3");
	set_variable_type				(3,SCALAR);
	set_variable_equation_type		(3,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(3, "c, n1, n2, n3, grad(u)");
    set_dependencies_gradient_term_RHS(3, "grad(n3)");


	// Variable 4
	set_variable_name				(4,"u");
	set_variable_type				(4,VECTOR);
	set_variable_equation_type		(4,TIME_INDEPENDENT);

    set_dependencies_value_term_RHS(4, "");
    set_dependencies_gradient_term_RHS(4, "c, n1, n2, n3, grad(u)");
    set_dependencies_value_term_LHS(4, "");
    set_dependencies_gradient_term_LHS(4, "grad(change(u))");

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

#define h1V (n1*n1*n1*(10.0-15.0*n1+6.0*n1*n1))
#define hn1V (30.0*(n1-1.0)*(n1-1.0)*n1*n1)
#define h2V (n2*n2*n2*(10.0-15.0*n2+6.0*n2*n2))
#define hn2V (30.0*(n2-1.0)*(n2-1.0)*n2*n2)
#define h3V (n3*n3*n3*(10.0-15.0*n3+6.0*n3*n3))
#define hn3V (30.0*(n3-1.0)*(n3-1.0)*n3*n3)

// #define h1V (3.0*n1*n1-2.0*n1*n1*n1)
// #define hn1V (6.0*n1-6.0*n1*n1)
// #define h2V (3.0*n2*n2-2.0*n2*n2*n2)
// #define hn2V (6.0*n2-6.0*n2*n2)
// #define h3V (3.0*n3*n3-2.0*n3*n3*n3)
// #define hn3V (6.0*n3-6.0*n3*n3)

// This double-well function can be used to tune the interfacial energy
#define fbarrierV (n1*n1-2.0*n1*n1*n1+n1*n1*n1*n1 + n2*n2-2.0*n2*n2*n2+n2*n2*n2*n2 + n3*n3-2.0*n3*n3*n3+n3*n3*n3*n3 + 5.0*(n1*n1*n2*n2 + n1*n1*n3*n3 + n2*n2*n3*n3) + 5.0*n1*n1*n2*n2*n3*n3)
#define fbarriern1V (2.0*n1-6.0*n1*n1+4.0*n1*n1*n1 + 10.0*n1*(n2*n2+n3*n3) + 10.0*n1*n2*n2*n3*n3)
#define fbarriern2V (2.0*n2-6.0*n2*n2+4.0*n2*n2*n2 + 10.0*n2*(n1*n1+n3*n3) + 10.0*n2*n1*n1*n3*n3)
#define fbarriern3V (2.0*n3-6.0*n3*n3+4.0*n3*n3*n3 + 10.0*n3*(n2*n2+n1*n1) + 10.0*n3*n2*n2*n1*n1)

// Residuals
#define rcV   (c)
#define rcxTemp ( cx + (c_alpha-c_beta)*(hn1V*n1x + hn2V*n2x + hn3V*n3x))
#define rcxV  (constV(-userInputs.dtValue)*McV*rcxTemp)

#define rn1V   (n1-constV(userInputs.dtValue*Mn1V)*gamma*( (fbV-faV)*hn1V - (c_beta-c_alpha)*facV*hn1V + W*fbarriern1V + nDependentMisfitAC1))
#define rn1xV  (constV(-userInputs.dtValue*Mn1V)*gamma*Knx1)

#define rn2V   (n2-constV(userInputs.dtValue*Mn2V)*gamma*( (fbV-faV)*hn2V - (c_beta-c_alpha)*facV*hn2V + W*fbarriern2V + nDependentMisfitAC2))
#define rn2xV  (constV(-userInputs.dtValue*Mn2V)*gamma*Knx2)

#define rn3V   (n3-constV(userInputs.dtValue*Mn3V)*gamma*( (fbV-faV)*hn3V - (c_beta-c_alpha)*facV*hn3V + W*fbarriern3V + nDependentMisfitAC3))
#define rn3xV  (constV(-userInputs.dtValue*Mn3V)*gamma*Knx3)


// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---

// The concentration and its derivatives
scalarvalueType c = variable_list.get_scalar_value(0);
scalargradType cx = variable_list.get_scalar_gradient(0);

// The first order parameter and its derivatives
scalarvalueType n1 = variable_list.get_scalar_value(1);
scalargradType n1x = variable_list.get_scalar_gradient(1);

// The second order parameter and its derivatives
scalarvalueType n2 = variable_list.get_scalar_value(2);
scalargradType n2x = variable_list.get_scalar_gradient(2);

// The third order parameter and its derivatives
scalarvalueType n3 = variable_list.get_scalar_value(3);
scalargradType n3x = variable_list.get_scalar_gradient(3);

// The derivative of the displacement vector
vectorgradType ux = variable_list.get_vector_gradient(4);

vectorgradType ruxV;

// -------------------------------------------------
// Nucleation expressions
// -------------------------------------------------
dealii::AlignedVector<dealii::VectorizedArray<double > > source_terms(3,constV(0.0));
dealii::VectorizedArray<double> gamma = constV(1.0);
seedNucleus(q_point_loc,source_terms,gamma);

if (this->currentTime < 50){
    gamma = constV(0.0);
}

// -------------------------------------------------

// Calculate c_alpha and c_beta from c
scalarvalueType sum_hpV = h1V+h2V+h3V;
scalarvalueType c_alpha = ((B2*c+0.5*(B1-A1)*sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV));
scalarvalueType c_beta  = ((A2*c+0.5*(A1-B1)*(1.0-sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV)));

// Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
// Note: this section can be optimized to reduce recalculations (note: that doesn't change the performance)
scalarvalueType cbn1V, cbn2V, cbn3V, cbcV;


cbcV = faccV/( (constV(1.0)-sum_hpV)*fbccV + (sum_hpV)*faccV );

cbn1V = hn1V * (c_alpha - c_beta) * cbcV;
cbn2V = hn2V * (c_alpha - c_beta) * cbcV;
cbn3V = hn3V * (c_alpha - c_beta) * cbcV;

// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts1, sfts2, sfts3;

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	// Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c_beta + b_p
	sfts1[i][j] = constV(sfts_linear1[i][j])*c_beta + constV(sfts_const1[i][j]);
	sfts2[i][j] = constV(sfts_linear2[i][j])*c_beta + constV(sfts_const2[i][j]);
	sfts3[i][j] = constV(sfts_linear3[i][j])*c_beta + constV(sfts_const3[i][j]);
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
computeStress<dim>(CIJ_Mg, E2, S);

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


variable_list.set_scalar_value_term_RHS(0,rcV);
variable_list.set_scalar_gradient_term_RHS(0,rcxV);

variable_list.set_scalar_value_term_RHS(1,rn1V + source_terms[0]);
variable_list.set_scalar_gradient_term_RHS(1,rn1xV);

variable_list.set_scalar_value_term_RHS(2,rn2V + source_terms[1]);
variable_list.set_scalar_gradient_term_RHS(2,rn2xV);

variable_list.set_scalar_value_term_RHS(3,rn3V + source_terms[2]);
variable_list.set_scalar_gradient_term_RHS(3,rn3xV);

}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are not
// explicit time-dependent equations. It takes "variable_list" as an input, which is
// a list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one proportional
// to the test function and one proportional to the gradient of the test function. The
// index for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

 // --- Getting the values and derivatives of the model variables ---

 // The concentration and its derivatives
 scalarvalueType c = variable_list.get_scalar_value(0);

 // The first order parameter and its derivatives
 scalarvalueType n1 = variable_list.get_scalar_value(1);

 // The second order parameter and its derivatives
 scalarvalueType n2 = variable_list.get_scalar_value(2);

 // The third order parameter and its derivatives
 scalarvalueType n3 = variable_list.get_scalar_value(3);

 // The derivative of the displacement vector
 vectorgradType ux = variable_list.get_vector_gradient(4);

 // Calculate c_alpha and c_beta from c
 scalarvalueType sum_hpV = h1V+h2V+h3V;
 scalarvalueType c_alpha = ((B2*c+0.5*(B1-A1)*sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV));
 scalarvalueType c_beta  = ((A2*c+0.5*(A1-B1)*(1.0-sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV)));


 // Calculate the stress-free transformation strain and its derivatives at the quadrature point
 dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts1, sfts2, sfts3;

 for (unsigned int i=0; i<dim; i++){
 for (unsigned int j=0; j<dim; j++){
 	// Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c_beta + b_p
 	sfts1[i][j] = constV(sfts_linear1[i][j])*c_beta + constV(sfts_const1[i][j]);
 	sfts2[i][j] = constV(sfts_linear2[i][j])*c_beta + constV(sfts_const2[i][j]);
 	sfts3[i][j] = constV(sfts_linear3[i][j])*c_beta + constV(sfts_const3[i][j]);
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
 computeStress<dim>(CIJ_Mg, E2, S);

 // Fill residual corresponding to mechanics
 // R=-C*(E-E0)
 vectorgradType ruxV;
 for (unsigned int i=0; i<dim; i++){
 for (unsigned int j=0; j<dim; j++){
 	  ruxV[i][j] = - S[i][j];
 }
 }

 variable_list.set_vector_gradient_term_RHS(4,ruxV);

}

// =================================================================================
// seedNucleus: a function particular to this app
// =================================================================================
template <int dim,int degree>
void customPDE<dim,degree>::seedNucleus(const dealii::Point<dim, dealii::VectorizedArray<double> > & q_point_loc,
	//std::vector<dealii::VectorizedArray<double> > & source_terms,
        dealii::AlignedVector<dealii::VectorizedArray<double> > & source_terms,
	dealii::VectorizedArray<double> & gamma) const {

		for (typename std::vector<nucleus<dim> >::const_iterator thisNucleus=this->nuclei.begin(); thisNucleus!=this->nuclei.end(); ++thisNucleus){

			if (thisNucleus->seededTime + thisNucleus->seedingTime > this->currentTime){

				// Calculate the weighted distance function to the order parameter freeze boundary (weighted_dist = 1.0 on that boundary)
				dealii::VectorizedArray<double> weighted_dist = this->weightedDistanceFromNucleusCenter(thisNucleus->center,
                    userInputs.get_nucleus_freeze_semiaxes(thisNucleus->orderParameterIndex),
                    userInputs.get_nucleus_rotation_matrix(thisNucleus->orderParameterIndex),
                    q_point_loc,
                    thisNucleus->orderParameterIndex);

				for (unsigned i=0; i<gamma.n_array_elements;i++){
					if (weighted_dist[i] <= 1.0){
						gamma[i] = 0.0;

						// Seed a nucleus if it was added to the list of nuclei this time step
						if (thisNucleus->seedingTimestep == this->currentIncrement){

                            // Set the rotation matrix for this nucleus (two possible per order parameter)
                            dealii::Tensor<2,dim,double> nucleus_rot_matrix;
                            if (thisNucleus->random_number < 0.5){
                                //std::cout << "a" << thisNucleus->orderParameterIndex << std::endl;
                                nucleus_rot_matrix = userInputs.get_nucleus_rotation_matrix(thisNucleus->orderParameterIndex);
                            }
                            else {
                                //std::cout << "b" << thisNucleus->orderParameterIndex << std::endl;
                                double degrees_to_rad = std::acos(0.0)/90.0;

                                dealii::Tensor<2,dim,double> Rz;
                                Rz[0][0] = std::cos(-49.4*degrees_to_rad);
                                Rz[1][0] = std::sin(-49.4*degrees_to_rad);
                                Rz[0][1] = -std::sin(-49.4*degrees_to_rad);
                                Rz[1][1] = std::cos(-49.4*degrees_to_rad);

                                if (dim == 3){
                                    Rz[2][2] = 1.0;
                                }
                                nucleus_rot_matrix = userInputs.get_nucleus_rotation_matrix(thisNucleus->orderParameterIndex)*Rz;
                            }

							// Find the weighted distance to the outer edge of the nucleus and use it to calculate the order parameter source term (r = 1.0 on that boundary)
							dealii::Point<dim,double> q_point_loc_element;
							for (unsigned int j=0; j<dim; j++){
								q_point_loc_element(j) = q_point_loc(j)[i];
							}
							double r = this->weightedDistanceFromNucleusCenter(thisNucleus->center,
                                userInputs.get_nucleus_semiaxes(thisNucleus->orderParameterIndex),
                                nucleus_rot_matrix,
                                q_point_loc_element,
                                thisNucleus->orderParameterIndex);

							double avg_semiaxis = 0.0;
							for (unsigned int j=0; j<dim; j++){
								avg_semiaxis += thisNucleus->semiaxes[j];
							}
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

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of that
// quadrature point is given by "q_point_loc". The function outputs two terms to
// variable_list -- one proportional to the test function and one proportional to the
// gradient of the test function -- for the left-hand-side of the equation. The index
// for each variable in this list corresponds to the index given at the top of this
// file. If there are multiple elliptic equations, conditional statements should be
// sed to ensure that the correct residual is being submitted. The index of the field
// being solved can be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

    vectorgradType ruxV;

    // Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the dealii "symmetrize" function
    dealii::Tensor<2, dim, dealii::VectorizedArray<double> > E = symmetrize(variable_list.get_change_in_vector_gradient(4));

    // Compute stress tensor (which is equal to the residual, Rux)
    computeStress<dim>(CIJ_Mg, E,ruxV);

    variable_list.set_vector_gradient_term_LHS(4,ruxV);

}
