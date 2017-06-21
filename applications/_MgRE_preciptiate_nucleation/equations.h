// List of variables and residual equations for the Precipitate Evolution example application

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".

// // Cahn-Hilliard mobility
// #define McV 1.0
//
// // Allen-Cahn mobilities
// #define Mn1V (100.0)
// #define Mn2V (100.0)
// #define Mn3V (100.0)
//
// // Gradient energy coefficients
// //double Kn1[3][3]={{0.01141*scaleFactor,0,0},{0,0.01426*scaleFactor,0},{0,0,0.004326*scaleFactor}}; // Scaled KKS B'''
// double Kn1[3][3]={{0.0318/1.5,0,0},{0,(6.7347e-4)/1.5,0},{0,0,0.0269/1.5}}; // Scaled KKS B' Mg-Y
// //double Kn1[3][3]={{0.01141*scaleFactor,0,0},{0,0.01141*scaleFactor,0},{0,0,0.01141*scaleFactor}}; // Isotropic interfacial energy
//
// double Kn2[3][3]={{0.008455/1.5,-0.01348/1.5,0},{-0.01348/1.5,0.02402/1.5,0},{0,0,0.0269/1.5}}; // Scaled KKS B'''
// double Kn3[3][3]={{0.008455/1.5,0.01348/1.5,0},{0.01348/1.5,0.02402/1.5,0},{0,0,0.0269/1.5}}; // Scaled KKS B'''
// //define energy barrier coefficient (used to tune the interfacial energy)
// #define W (0.1288*1.5)
//
// // Define Mechanical properties
// #define n_dependent_stiffness false
//
// // Stress-free transformation strains
// // Linear fits for the stress-free transformation strains in for sfts = sfts_linear * c + sfts_const
//
// // B''' (gen 2)
// //double sfts_linear1[3][3] = {{-0.34358,0,0},{0,0.68568,0},{0,0,0.19308}};
// //double sfts_const1[3][3] = {{0.14978,0,0},{0,-0.10254,0},{0,0,-0.034049}};
//
// // B''' (gen 4)
// //double sfts_linear1[3][3] = {{-0.32067,0,0},{0,0.66232,0},{0,0,0.19462}};
// //double sfts_const1[3][3] = {{0.14698,0,0},{0,-0.09877,0},{0,0,-0.034899}};
// //
// //double sfts_linear2[3][3] = {{0.4166,0.4256,0},{0.4256,-0.07495,0},{0,0,0.19462}};
// //double sfts_const2[3][3] = {{-0.03733,-0.1064,0},{-0.1064,0.08554,0},{0,0,-0.034899}};
// //
// //double sfts_linear3[3][3] = {{0.4166,-0.4256,0},{-0.4256,-0.07495,0},{0,0,0.19462}};
// //double sfts_const3[3][3] = {{-0.03733,0.1064,0},{0.1064,0.08554,0},{0,0,-0.034899}};
//
// // B's (Mg-Y)
// double sfts_linear1[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
// double sfts_const1[3][3] = {{(0.013972+0.13383*0.125),0,0},{0,(-0.0057545+0.21273*0.125),0},{0,0,-0.002691+0.014832*0.125}};
// double sfts_linear2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
// double sfts_const2[3][3] = {{0.02330,-0.004271,0},{-0.004271,0.0282,0},{0,0,-0.0008370}};
// double sfts_linear3[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
// double sfts_const3[3][3] = {{0.02330,0.004271,0},{0.004271,0.0282,0},{0,0,-0.0008370}};
//
//
// // Zero misfit
// //double sfts_linear1[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
// //double sfts_const1[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
// //double sfts_linear2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
// //double sfts_const2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
// //double sfts_linear3[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
// //double sfts_const3[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
//
// ////define free energy expressions (Mg-Nd data from CASM) (B''', gen 4)
// //double A2 = 100.56;
// //double A1 = -1.727;
// //double A0 = 0.0001138;
// ////double B2 = 2.4929;
// ////double B1 = -2.2810;
// ////double B0 = 0.039048;
// //
// //// B''' fits without B''
// //double B2 = 5.2408;
// //double B1 = -2.9679;
// //double B0 = 0.081978;
//
// // Stand-in Mg-Y free energies
// //define free energy expressions (Mg-Nd data from CASM) (B''', gen 4)
// double A2 = 7.5096;
// double A1 = -0.0078096;
// double A0 = 2.0305e-6;
// double B2 = 2.8557;
// double B1 = -0.71393;
// double B0 = 0.04462;

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
// Set the nucleation parameters
// =================================================================================

#define interface_coeff 0.5

#define epsil 1.0e-7

// Constants k1 and k2 for nucleation rate in the bulk
#define k1 (8.39e-11*1000000.0)
#define k2 4.01e-4
#define tau (7523.0*0.0)



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
template <int dim,int degree>
void customPDE<dim,degree>::residualRHS(const std::vector<modelVariable<dim> > & modelVariablesList,
												std::vector<modelResidual<dim> > & modelResidualsList,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// The concentration and its derivatives (names here should match those in the macros above)
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

// The first order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n1 = modelVariablesList[1].scalarValue;
scalargradType n1x = modelVariablesList[1].scalarGrad;

// The second order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n2 = modelVariablesList[2].scalarValue;
scalargradType n2x = modelVariablesList[2].scalarGrad;

// The third order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n3 = modelVariablesList[3].scalarValue;
scalargradType n3x = modelVariablesList[3].scalarGrad;

// The derivative of the displacement vector (names here should match those in the macros above)
vectorgradType ux = modelVariablesList[4].vectorGrad;
vectorgradType ruxV;

vectorhessType uxx;

bool c_dependent_misfit = false;
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		if (std::abs(sfts_linear1[i][j])>1.0e-12){
			c_dependent_misfit = true;
		}
	}
}

if (c_dependent_misfit == true){
	uxx = modelVariablesList[4].vectorHess;
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
		  CIJ_combined[i][j] = this->userInputs.CIJ_list[0][i][j]*(constV(1.0)-sum_hpV) + this->userInputs.material_moduli.CIJ_list[1][i][j]*(h1V) + this->userInputs.material_moduli.CIJ_list[2][i][j]*(h2V) + this->userInputs.material_moduli.CIJ_list[3][i][j]*(h3V);
	  }
}
computeStress<dim>(CIJ_combined, E2, S);
}
else{
computeStress<dim>(this->userInputs.material_moduli.CIJ_list[0], E2, S);
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
	computeStress<dim>(this->userInputs.material_moduli.CIJ_list[1]-this->userInputs.material_moduli.CIJ_list[0], E2, S2);
	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			heterMechAC1 += S2[i][j]*E2[i][j];
		}
	}
	heterMechAC1 *= 0.5*hn1V;

	computeStress<dim>(this->userInputs.material_moduli.CIJ_list[2]-this->userInputs.material_moduli.CIJ_list[0], E2, S2);
	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			heterMechAC2 += S2[i][j]*E2[i][j];
		}
	}
	heterMechAC2 *= 0.5*hn2V;

	computeStress<dim>(this->userInputs.material_moduli.CIJ_list[3]-this->userInputs.material_moduli.CIJ_list[0], E2, S2);
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
		computeStress<dim>(this->userInputs.material_moduli.CIJ_list[0], E3, S3);
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
					computeStress<dim>(this->userInputs.material_moduli.CIJ_list[1]-this->userInputs.material_moduli.CIJ_list[0], E2, S2);
					grad_mu_el[k]+= - S2[i][j] * (cbcV*E4[i][j]*hn1V*n1x[k]);
					computeStress<dim>(this->userInputs.material_moduli.CIJ_list[2]-this->userInputs.material_moduli.CIJ_list[0], E2, S2);
					grad_mu_el[k]+= - S2[i][j] * (cbcV*E4[i][j]*hn2V*n2x[k]);
					computeStress<dim>(this->userInputs.material_moduli.CIJ_list[3]-this->userInputs.material_moduli.CIJ_list[0], E2, S2);
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

modelResidualsList[0].scalarValueResidual = rcV;
modelResidualsList[0].scalarGradResidual = rcxV;

modelResidualsList[1].scalarValueResidual = rn1V + source_terms[0];
modelResidualsList[1].scalarGradResidual = rn1xV;
modelResidualsList[2].scalarValueResidual = rn2V + source_terms[1];
modelResidualsList[2].scalarGradResidual = rn2xV;
modelResidualsList[3].scalarValueResidual = rn3V + source_terms[2];
modelResidualsList[3].scalarGradResidual = rn3xV;

modelResidualsList[4].vectorGradResidual = ruxV;

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

				// Calculate the weighted distance function to the order parameter freeze boundary (weighted_dist = 1.0 on that boundary)
				dealii::VectorizedArray<double> weighted_dist = constV(0.0);
				for (unsigned int i=0; i<dim; i++){
					dealii::VectorizedArray<double> temp = (thisNucleus->center(i) - q_point_loc(i));

					if (userInputs.BC_list[1].var_BC_type[2*i]==PERIODIC){
						for (unsigned j=0; j<gamma.n_array_elements;j++)
						temp[j] -= round(temp[j]/userInputs.domain_size[i])*userInputs.domain_size[i];
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
							double r = 0.0;
							double avg_semiaxis = 0.0;
							for (unsigned int j=0; j<dim; j++){
								double temp = (thisNucleus->center(j) - q_point_loc(j)[i]);
								if (userInputs.BC_list[1].var_BC_type[2*i]==PERIODIC){
									double domsize_j = userInputs.domain_size[j];
									temp=temp-round(temp/domsize_j)*domsize_j;
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
template <int dim,int degree>
void customPDE<dim,degree>::residualLHS(const std::vector<modelVariable<dim> > & modelVariablesList,
		modelResidual<dim> & modelRes,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//n1
scalarvalueType n1 = modelVariablesList[0].scalarValue;

//n2
scalarvalueType n2 = modelVariablesList[1].scalarValue;

//n3
scalarvalueType n3 = modelVariablesList[2].scalarValue;

// Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the dealii "symmetrize" function
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > E = symmetrize(modelVariablesList[3].vectorGrad);

// Compute stress tensor (which is equal to the residual, Rux)
if (n_dependent_stiffness == true){
	dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > CIJ_combined;
	CIJ_combined = this->userInputs.material_moduli.CIJ_list[0]*(constV(1.0)-h1V-h2V-h3V) + this->userInputs.material_moduli.CIJ_list[1]*(h1V) + this->userInputs.material_moduli.CIJ_list[2]*(h2V) + this->userInputs.material_moduli.CIJ_list[3]*(h3V);

	computeStress<dim>(CIJ_combined, E, modelRes.vectorGradResidual);
}
else{
	computeStress<dim>(this->userInputs.material_moduli.CIJ_list[0], E, modelRes.vectorGradResidual);
}

}

// =================================================================================
// energyDensity (needed only if calcEnergy == true)
// =================================================================================
// This function integrates the free energy density across the computational domain.
// It takes "modelVariablesList" as an input, which is a list of the value and
// derivatives of each of the variables at a specific quadrature point. It also
// takes the mapped quadrature weight, "JxW_value", as an input. The (x,y,z) location
// of the quadrature point is given by "q_point_loc". The weighted value of the
// energy density is added to "energy" variable and the components of the energy
// density are added to the "energy_components" variable (index 0: chemical energy,
// index 1: gradient energy, index 2: elastic energy).
template <int dim,int degree>
void customPDE<dim,degree>::energyDensity(const std::vector<modelVariable<dim> > & modelVariablesList,
											const dealii::VectorizedArray<double> & JxW_value,
											dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {

scalarvalueType total_energy_density = constV(0.0);

//c
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

//n1
scalarvalueType n1 = modelVariablesList[1].scalarValue;
scalargradType n1x = modelVariablesList[1].scalarGrad;

// The second order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n2 = modelVariablesList[2].scalarValue;
scalargradType n2x = modelVariablesList[2].scalarGrad;

// The third order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n3 = modelVariablesList[3].scalarValue;
scalargradType n3x = modelVariablesList[3].scalarGrad;

//u
vectorgradType ux = modelVariablesList[4].vectorGrad;

scalarvalueType sum_hpV = h1V+h2V+h3V;
scalarvalueType c_alpha = ((B2*c+0.5*(B1-A1)*sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV));
scalarvalueType c_beta  = ((A2*c+0.5*(A1-B1)*(1.0-sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV)));

scalarvalueType f_chem = (constV(1.0)-(h1V+h2V+h3V))*faV + (h1V+h2V+h3V)*fbV + W*(fbarrierV);

scalarvalueType f_grad = constV(0.0);

for (int i=0; i<dim; i++){
  for (int j=0; j<dim; j++){
	  f_grad += constV(0.5*Kn1[i][j])*n1x[i]*n1x[j] + constV(0.5*Kn2[i][j])*n2x[i]*n2x[j] + constV(0.5*Kn3[i][j])*n3x[i]*n3x[j];
  }
}


// Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
// Note: this section can be optimized to reduce recalculations
scalarvalueType cbn1V, cbn2V, cbn3V, cbcV, cbcn1V, cbcn2V, cbcn3V;

cbcV = faccV/( (constV(1.0)-sum_hpV)*fbccV + (sum_hpV)*faccV );

cbn1V = hn1V * (c_alpha - c_beta) * cbcV;
cbn2V = hn2V * (c_alpha - c_beta) * cbcV;
cbn3V = hn3V * (c_alpha - c_beta) * cbcV;

cbcn1V = (faccV * (fbccV-faccV) * hn1V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic
cbcn2V = (faccV * (fbccV-faccV) * hn2V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic
cbcn3V = (faccV * (fbccV-faccV) * hn3V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic


// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts1, sfts1c, sfts1cc, sfts2, sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc, sfts1n, sfts1cn, sfts2n, sfts2cn, sfts3n, sfts3cn;

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	// Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c_beta + b_p
	sfts1[i][j] = constV(sfts_linear1[i][j])*c_beta + constV(sfts_const1[i][j]);
	sfts1c[i][j] = constV(sfts_linear1[i][j]) * cbcV;
	sfts1cc[i][j] = constV(0.0);
	sfts1n[i][j] = constV(sfts_linear1[i][j]) * cbn1V;
	sfts1cn[i][j] = constV(sfts_linear1[i][j]) * cbcn1V;

	sfts2[i][j] = constV(sfts_linear2[i][j])*c_beta + constV(sfts_const2[i][j]);
	sfts2c[i][j] = constV(sfts_linear2[i][j]) * cbcV;
	sfts2cc[i][j] = constV(0.0);
	sfts2n[i][j] = constV(sfts_linear2[i][j]) * cbn2V;
	sfts2cn[i][j] = constV(sfts_linear2[i][j]) * cbcn2V;

	sfts3[i][j] = constV(sfts_linear3[i][j])*c_beta + constV(sfts_const3[i][j]);
	sfts3c[i][j] = constV(sfts_linear3[i][j]) * cbcV;
	sfts3cc[i][j] = constV(0.0);
	sfts3n[i][j] = constV(sfts_linear3[i][j]) * cbn3V;
	sfts3cn[i][j] = constV(sfts_linear3[i][j]) * cbcn2V;

}
}


//compute E2=(E-E0)
dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V + sfts2[i][j]*h2V + sfts3[i][j]*h3V);

  }
}

//compute stress
//S=C*(E-E0)
dealii::VectorizedArray<double> CIJ_combined[2*dim-1+dim/3][2*dim-1+dim/3];

if (n_dependent_stiffness == true){
  for (unsigned int i=0; i<2*dim-1+dim/3; i++){
	  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
		  CIJ_combined[i][j] = this->userInputs.material_moduli.CIJ_list[0][i][j]*(constV(1.0)-h1V-h2V-h3V) + this->userInputs.material_moduli.CIJ_list[1][i][j]*(h1V) + this->userInputs.material_moduli.CIJ_list[2][i][j]*(h2V) + this->userInputs.material_moduli.CIJ_list[3][i][j]*(h3V);
	  }
  }
  computeStress<dim>(CIJ_combined, E2, S);
}
else{
  computeStress<dim>(this->userInputs.material_moduli.CIJ_list[0], E2, S);
}

scalarvalueType f_el = constV(0.0);

for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  f_el += constV(0.5) * S[i][j]*E2[i][j];
  }
}

total_energy_density = f_chem + f_grad + f_el;

// Loop to step through each element of the vectorized arrays. Working with deal.ii
// developers to see if there is a more elegant way to do this.
this->assembler_lock.acquire ();
for (unsigned i=0; i<c.n_array_elements;i++){
  if (c[i] > 1.0e-10){
	  this->energy+=total_energy_density[i]*JxW_value[i];
	  this->energy_components[0]+= f_chem[i]*JxW_value[i];
	  this->energy_components[1]+= f_grad[i]*JxW_value[i];
	  this->energy_components[2]+= f_el[i]*JxW_value[i];
  }
}
this->assembler_lock.release ();
}
