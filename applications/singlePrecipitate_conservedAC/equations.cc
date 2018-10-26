// List of variables and residual equations for the Precipitate Evolution example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"c");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "c, n1, grad(u)");
    set_dependencies_gradient_term_RHS(0, "");

	// Variable 1
	set_variable_name				(1,"n1");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(1, "c, n1, grad(u)");
    set_dependencies_gradient_term_RHS(1, "grad(n1)");

	// Variable 2
	set_variable_name				(2,"u");
	set_variable_type				(2,VECTOR);
	set_variable_equation_type		(2,TIME_INDEPENDENT);

    set_dependencies_value_term_RHS(2, "");
    set_dependencies_gradient_term_RHS(2, "c, n1, grad(u)");
    set_dependencies_value_term_LHS(2, "");
    set_dependencies_gradient_term_LHS(2, "n1, grad(change(u))");

}

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".

// Calculate c_alpha and c_beta from c
#define c_alpha ((B2*c+0.5*(B1-A1)*h1V)/(A2*h1V+B2*(1.0-h1V)))
#define c_beta ((A2*c+0.5*(A1-B1)*(1.0-h1V))/(A2*h1V+B2*(1.0-h1V)))

#define faV (A2*c_alpha*c_alpha + A1*c_alpha + A0)
#define facV (2.0*A2*c_alpha +A1)
#define faccV (2.0*A2)
#define fbV (B2*c_beta*c_beta + B1*c_beta + B0)
#define fbcV (2.0*B2*c_beta +B1)
#define fbccV (2.0*B2)

#define h1V (3.0*n1*n1-2.0*n1*n1*n1)
#define hn1V (6.0*n1-6.0*n1*n1)

// This double-well function can be used to tune the interfacial energy
#define fbarrierV (n1*n1-2.0*n1*n1*n1+n1*n1*n1*n1)
#define fbarriernV (2.0*n1-6.0*n1*n1+4.0*n1*n1*n1)

// Residuals
#define rcV   (c-constV(userInputs.dtValue*dt_modifier*McV)*mu_c)

#define rn1V   (n1-constV(userInputs.dtValue*dt_modifier*Mn1V)*( (fbV-faV)*hn1V - (c_beta-c_alpha)*facV*hn1V + W*fbarriernV + nDependentMisfitAC1 + heterMechAC1))
#define rn1xV  (constV(-userInputs.dtValue*dt_modifier*Mn1V)*Knx1)


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

// The concentration and its derivatives (names here should match those in the macros above)
scalarvalueType c = variable_list.get_scalar_value(0);

// The first order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n1 = variable_list.get_scalar_value(1);
scalargradType n1x = variable_list.get_scalar_gradient(1);

// The derivative of the displacement vector (names here should match those in the macros above)
vectorgradType ux = variable_list.get_vector_gradient(2);


// Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
scalarvalueType cbnV, cbcV, cacV;

cbcV = faccV/( (constV(1.0)-h1V)*fbccV + h1V*faccV );
cacV = fbccV/( (constV(1.0)-h1V)*fbccV + h1V*faccV );
cbnV = hn1V * (c_alpha - c_beta) * cbcV;

// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts1;

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	// Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c_beta + b_p
	sfts1[i][j] = constV(sfts_const1[i][j]);
}
}

//compute E2=(E-E0)
dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V );
}
}

//compute stress
//S=C*(E-E0)
// Compute stress tensor (which is equal to the residual, Rux)
dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

if (n_dependent_stiffness == true){
for (unsigned int i=0; i<2*dim-1+dim/3; i++){
	  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
		  CIJ_combined[i][j] = CIJ_Mg[i][j]*(constV(1.0)-h1V) + CIJ_Beta[i][j]*h1V;
	  }
}
computeStress<dim>(CIJ_combined, E2, S);
}
else{
computeStress<dim>(CIJ_Mg, E2, S);
}


// Compute one of the stress terms in the order parameter chemical potential, nDependentMisfitACp = -C*(E-E0)*(E0_n)
dealii::VectorizedArray<double> nDependentMisfitAC1=constV(0.0);

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  nDependentMisfitAC1+=-S[i][j]*(sfts1[i][j]*hn1V);
}
}

// Compute the other stress term in the order parameter chemical potential, heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
dealii::VectorizedArray<double> heterMechAC1=constV(0.0);
dealii::VectorizedArray<double> S2[dim][dim];

if (n_dependent_stiffness == true){
	computeStress<dim>(CIJ_Beta-CIJ_Mg, E2, S2);

	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			heterMechAC1 += S2[i][j]*E2[i][j];
		}
	}
	heterMechAC1 = 0.5*hn1V*heterMechAC1;
}

//compute K*nx
scalargradType Knx1;
for (unsigned int a=0; a<dim; a++) {
Knx1[a]=0.0;
for (unsigned int b=0; b<dim; b++){
	  Knx1[a]+=constV(Kn1[a][b])*n1x[b];
}
}

// Calculate the chemical potential for the concentration
scalarvalueType mu_c = constV(0.0);
mu_c += facV*cacV * (1.0-h1V) + fbcV*cbcV * h1V;


variable_list.set_scalar_value_term_RHS(0,rcV);

variable_list.set_scalar_value_term_RHS(1,rn1V);
variable_list.set_scalar_gradient_term_RHS(1,rn1xV);

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

 // The concentration and its derivatives (names here should match those in the macros above)
 scalarvalueType c = variable_list.get_scalar_value(0);

 // The first order parameter and its derivatives (names here should match those in the macros above)
 scalarvalueType n1 = variable_list.get_scalar_value(1);

 // The derivative of the displacement vector (names here should match those in the macros above)
 vectorgradType ux = variable_list.get_vector_gradient(2);

 vectorgradType ruxV;

 // Calculate the stress-free transformation strain and its derivatives at the quadrature point
 dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts1;

 for (unsigned int i=0; i<dim; i++){
 for (unsigned int j=0; j<dim; j++){
 	// Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c_beta + b_p
 	sfts1[i][j] = constV(sfts_const1[i][j]);
 }
 }

 //compute E2=(E-E0)
 dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

 for (unsigned int i=0; i<dim; i++){
 for (unsigned int j=0; j<dim; j++){
 	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V );
 }
 }

 //compute stress
 //S=C*(E-E0)
 // Compute stress tensor (which is equal to the residual, Rux)
 dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

 if (n_dependent_stiffness == true){
 for (unsigned int i=0; i<2*dim-1+dim/3; i++){
 	  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
 		  CIJ_combined[i][j] = CIJ_Mg[i][j]*(constV(1.0)-h1V) + CIJ_Beta[i][j]*h1V;
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

 variable_list.set_vector_gradient_term_RHS(2,ruxV);

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

	//n1
	scalarvalueType n1 = variable_list.get_scalar_value(1);

	vectorgradType ruxV;

	// Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the dealii "symmetrize" function
	dealii::Tensor<2, dim, dealii::VectorizedArray<double> > E;
	E = symmetrize(variable_list.get_change_in_vector_gradient(2));

	// Compute stress tensor (which is equal to the residual, Rux)
	if (n_dependent_stiffness == true){
		dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > CIJ_combined;
		CIJ_combined = CIJ_Mg*(constV(1.0)-h1V);
		CIJ_combined += CIJ_Beta*(h1V);

		computeStress<dim>(CIJ_combined, E, ruxV);
	}
	else{
		computeStress<dim>(CIJ_Mg, E, ruxV);
	}

	variable_list.set_vector_gradient_term_LHS(2,ruxV);

}
