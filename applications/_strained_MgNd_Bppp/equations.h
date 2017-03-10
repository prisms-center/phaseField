// List of variables and residual equations for the Precipitate Evolution example application

// =================================================================================
// Define the variables in the model
// =================================================================================
// The number of variables
#define num_var 3

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define variable_name {"c", "n1", "u"}
#define variable_type {"SCALAR","SCALAR","VECTOR"}
#define variable_eq_type {"PARABOLIC","PARABOLIC","ELLIPTIC"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define need_val {true, true, false}
#define need_grad {true, true, true}
#define need_hess {false, false, true}

// Flags for whether the residual equation has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define need_val_residual {true, true, false}
#define need_grad_residual {true, true, true}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqn
// for the left-hand-side of the iterative solver for elliptic equations
#define need_val_LHS {false, true, false}
#define need_grad_LHS {false, false, true}
#define need_hess_LHS {false, false, false, false}

// Flags for whether the residual equation for the left-hand-side of the iterative
// solver for elliptic equations has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define need_val_residual_LHS {false, false, false}
#define need_grad_residual_LHS {false, false, true}

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".

// Cahn-Hilliard mobility
#define McV 1.0

// Allen-Cahn mobilities
#define Mn1V (300.0/scaleFactor)

// Gradient energy coefficients
double Kn1[3][3]={{0.01141*scaleFactor,0,0},{0,0.01426*scaleFactor,0},{0,0,0.004326*scaleFactor}}; // Scaled KKS B'''
//double Kn1[3][3]={{0.01141*scaleFactor,0,0},{0,0.01141*scaleFactor,0},{0,0,0.01141*scaleFactor}}; // Isotropic interfacial energy

//define energy barrier coefficient (used to tune the interfacial energy)
#define W (0.8272/scaleFactor)

// Define Mechanical properties
#define n_dependent_stiffness true
// Mechanical symmetry of the material and stiffness parameters
#if problemDIM==1
	// Used throughout system if n_dependent_stiffness == false, used in n=0 phase if n_dependent_stiffness == true
	#define MaterialModelV ISOTROPIC
	#define MaterialConstantsV {22.5,0.3}
	// Used in n=1 phase if n_dependent_stiffness == true
	#define MaterialModelBetaV ISOTROPIC
	#define MaterialConstantsBetaV {22.5,0.3}

#elif problemDIM==2
	// If n_dependent_stiffness == false the first entry is used for all phases
	// 2D order of constants ANISOTROPIC - 6 constants [C11 C22 C33 C12 C13 C23]
	#define MaterialModels {{"ANISOTROPIC"},{"ANISOTROPIC"}}
	#define MaterialConstants {{31.3,31.3,6.65,13.0,0.0,0.0},{23.35,30.25,36.35,15.35,0.0,0.0}} //scaled by E* = 2e9 J/m^3

#elif problemDIM==3
// If n_dependent_stiffness == false the first entry is used for all phases
	#define MaterialModels {{"ANISOTROPIC"},{"ANISOTROPIC"}}
	// 3D order of constants ANISOTROPIC - 21 constants [11, 22, 33, 44, 55, 66, 12, 13, 14, 15, 16, 23, 24, 25, 26, 34, 35, 36, 45, 46, 56]
	//#define MaterialConstantsV {62.6,62.6,64.9,13.3,13.3,18.3,26.0,20.9,0,0,0,20.9,0,0,0,0,0,0,0,0,0} //these are in GPa-need to be non-dimensionalized
	#define MaterialConstants {{31.3,31.3,32.45,6.65,6.65,9.15,13.0,10.45,0,0,0,10.45,0,0,0,0,0,0,0,0,0},{23.35,30.25,36.35,8.2,16.7,14.45,15.35,14.35,0,0,0,7.25,0,0,0,0,0,0,0,0,0}} //scaled by E* = 2e9 J/m^3
#endif


// Stress-free transformation strains
// Linear fits for the stress-free transformation strains in for sfts = sfts_linear * c + sfts_const

// B''' (gen 2)
//double sfts_linear1[3][3] = {{-0.34358,0,0},{0,0.68568,0},{0,0,0.19308}};
//double sfts_const1[3][3] = {{0.14978,0,0},{0,-0.10254,0},{0,0,-0.034049}};

// B''' (gen 4)
double sfts_linear1[3][3] = {{-0.32067,0,0},{0,0.66232,0},{0,0,0.19462}};
double sfts_const1[3][3] = {{0.14698,0,0},{0,-0.09877,0},{0,0,-0.034899}};

// B'
//double sfts_linear1[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
//double sfts_const1[3][3] = {{-0.32067*0.125+0.14698,0,0},{0,0.66232*0.125-0.09877,0},{0,0,0.19462*0.125-0.034899}};

// Zero misfit
//double sfts_linear1[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
//double sfts_const1[3][3] = {{0,0,0},{0,0,0},{0,0,0}};


// Calculate c_alpha and c_beta from c
#define c_alpha ((B2*c+0.5*(B1-A1)*h1V)/(A2*h1V+B2*(1.0-h1V)))
#define c_beta ((A2*c+0.5*(A1-B1)*(1.0-h1V))/(A2*h1V+B2*(1.0-h1V)))

//define free energy expressions (Mg-Nd data from CASM) (B''', gen 4)
double A2 = 100.56;
double A1 = -1.727;
double A0 = 0.0001138;
//double B2 = 2.4929;
//double B1 = -2.2810;
//double B0 = 0.039048;

// B''' fits without B''
double B2 = 5.2408;
double B1 = -2.9679;
double B0 = 0.081978;

#define faV (A2*c_alpha*c_alpha + A1*c_alpha + A0)
#define facV (2.0*A2*c_alpha +A1)
#define faccV (2.0*A2)
#define fbV (B2*c_beta*c_beta + B1*c_beta + B0)
#define fbcV (2.0*B2*c_beta +B1)
#define fbccV (2.0*B2)

#define h1V (3.0*n1*n1-2.0*n1*n1*n1)
#define hn1V (6.0*n1-6.0*n1*n1)

#define h1strainV h1V
#define hn1strainV hn1V

// This double-well function can be used to tune the interfacial energy
#define fbarrierV (n1*n1-2.0*n1*n1*n1+n1*n1*n1*n1)
#define fbarriernV (2.0*n1-6.0*n1*n1+4.0*n1*n1*n1)

// Residuals
#define rcV   (c)
#define rcxTemp ( cx + n1x*(c_alpha-c_beta)*hn1V + grad_mu_el * (h1V*faccV+(constV(1.0)-h1V)*fbccV)/constV(faccV*fbccV) )
#define rcxV  (constV(-timeStep)*McV*rcxTemp)

#define rn1V   (n1-constV(timeStep*Mn1V)*( (fbV-faV)*hn1V - (c_beta-c_alpha)*facV*hn1V + W*fbarriernV + nDependentMisfitAC1 + heterMechAC1))
#define rn1xV  (constV(-timeStep*Mn1V)*Knx1)


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
void customPDE<dim,degree>::residualRHS(const std::vector<modelVariable<dim>> & modelVariablesList,
												std::vector<modelResidual<dim>> & modelResidualsList,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// The concentration and its derivatives (names here should match those in the macros above)
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

// The first order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n1 = modelVariablesList[1].scalarValue;
scalargradType n1x = modelVariablesList[1].scalarGrad;

// The derivative of the displacement vector (names here should match those in the macros above)
vectorgradType ux = modelVariablesList[2].vectorGrad;
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
	uxx = modelVariablesList[2].vectorHess;
}

// Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
scalarvalueType cbnV, cbcV, cbcnV;

cbcV = faccV/( (constV(1.0)-h1V)*fbccV + h1V*faccV );
cbnV = hn1V * (c_alpha - c_beta) * cbcV;
cbcnV = (faccV * (fbccV-faccV) * hn1V)/( ((1.0-h1V)*fbccV + h1V*faccV)*((1.0-h1V)*fbccV + h1V*faccV) );  // Note: this is only true if faV and fbV are quadratic

// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double> > sfts1, sfts1c, sfts1cc, sfts1n, sfts1cn;

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){

	// New SFTS formulation
	// Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c_beta + b_p
	sfts1[i][j] = constV(sfts_linear1[i][j])*c_beta + constV(sfts_const1[i][j]);
	sfts1c[i][j] = constV(sfts_linear1[i][j]) * cbcV;
	sfts1cc[i][j] = constV(0.0);
	sfts1n[i][j] = constV(sfts_linear1[i][j]) * cbnV;
	sfts1cn[i][j] = constV(sfts_linear1[i][j]) * cbcnV;
}
}

//compute E2=(E-E0)
dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1strainV );

}
}

//compute stress
//S=C*(E-E0)
// Compute stress tensor (which is equal to the residual, Rux)
dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

if (n_dependent_stiffness == true){
for (unsigned int i=0; i<2*dim-1+dim/3; i++){
	  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
		  CIJ_combined[i][j] = this->userInputs.CIJ_list[0][i][j]*(constV(1.0)-h1strainV) + this->userInputs.CIJ_list[1][i][j]*h1strainV;
	  }
}
computeStress<dim>(CIJ_combined, E2, S);
}
else{
computeStress<dim>(this->userInputs.CIJ_list[0], E2, S);
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

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  nDependentMisfitAC1+=-S[i][j]*(sfts1n[i][j]*h1V + sfts1[i][j]*hn1V);
}
}


// Compute the other stress term in the order parameter chemical potential, heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
dealii::VectorizedArray<double> heterMechAC1=constV(0.0);
dealii::VectorizedArray<double> S2[dim][dim];

if (n_dependent_stiffness == true){
//computeStress<dim>(CIJ_diff, E2, S2);
	computeStress<dim>(this->userInputs.CIJ_list[1]-this->userInputs.CIJ_list[0], E2, S2);
for (unsigned int i=0; i<dim; i++){
	  for (unsigned int j=0; j<dim; j++){
		  heterMechAC1 += S2[i][j]*E2[i][j];
	  }
}

heterMechAC1 = 0.5*hn1strainV*heterMechAC1;
}

// compute the stress term in the gradient of the concentration chemical potential, grad_mu_el = -[C*(E-E0)*E0c]x, must be a vector with length dim
scalargradType grad_mu_el;

if (c_dependent_misfit == true){
	dealii::VectorizedArray<double> E3[dim][dim], S3[dim][dim];

	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			E3[i][j] =  -( sfts1c[i][j]*h1strainV );
		}
	}

	if (n_dependent_stiffness == true){
		computeStress<dim>(CIJ_combined, E3, S3);
	}
	else{
		computeStress<dim>(this->userInputs.CIJ_list[0], E3, S3);
	}

	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			for (unsigned int k=0; k<dim; k++){
				grad_mu_el[k] += S3[i][j] * (constV(0.5)*(uxx[i][j][k]+uxx[j][i][k]) + E3[i][j]*cx[k]
														  - ( (sfts1[i][j]*hn1strainV + sfts1n[i][j]*h1V) *n1x[k]) );

				grad_mu_el[k]+= - S[i][j] * ( (sfts1c[i][j]*hn1strainV + h1V*sfts1cn[i][j])*n1x[k] + ( sfts1cc[i][j]*h1strainV )*cx[k]);

				if (n_dependent_stiffness == true){
					grad_mu_el[k]+= S2[i][j] * (hn1strainV*n1x[k]) * (E3[i][j]);

				}
			}
		}
	}
}


//compute K*nx
scalargradType Knx1;
for (unsigned int a=0; a<dim; a++) {
Knx1[a]=0.0;
for (unsigned int b=0; b<dim; b++){
	  Knx1[a]+=constV(Kn1[a][b])*n1x[b];
}
}

modelResidualsList[0].scalarValueResidual = rcV;
modelResidualsList[0].scalarGradResidual = rcxV;

modelResidualsList[1].scalarValueResidual = rn1V;
modelResidualsList[1].scalarGradResidual = rn1xV;

modelResidualsList[2].vectorGradResidual = ruxV;

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
void customPDE<dim,degree>::residualLHS(const std::vector<modelVariable<dim>> & modelVariablesList,
		modelResidual<dim> & modelRes,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//n1
scalarvalueType n1 = modelVariablesList[0].scalarValue;

//u
vectorgradType ux = modelVariablesList[1].vectorGrad;
vectorgradType ruxV;

// Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the dealii "symmetrize" function
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > E;
E = symmetrize(ux);

// Compute stress tensor (which is equal to the residual, Rux)
if (n_dependent_stiffness == true){
	dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > CIJ_combined;
	CIJ_combined = this->userInputs.CIJ_list[0]*(constV(1.0)-h1strainV);
	CIJ_combined += this->userInputs.CIJ_list[1]*(h1strainV);

	computeStress<dim>(CIJ_combined, E, ruxV);
}
else{
	computeStress<dim>(this->userInputs.CIJ_list[0], E, ruxV);
}

modelRes.vectorGradResidual = ruxV;

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
template <int dim, int degree>
void customPDE<dim,degree>::energyDensity(const std::vector<modelVariable<dim>> & modelVarList,
											const dealii::VectorizedArray<double> & JxW_value,
											dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {

scalarvalueType total_energy_density = constV(0.0);

//c
scalarvalueType c = modelVarList[0].scalarValue;
scalargradType cx = modelVarList[0].scalarGrad;

//n1
scalarvalueType n1 = modelVarList[1].scalarValue;
scalargradType n1x = modelVarList[1].scalarGrad;

//u
vectorgradType ux = modelVarList[2].vectorGrad;

scalarvalueType f_chem = (constV(1.0)-(h1V))*faV + (h1V)*fbV;

scalarvalueType f_grad = constV(0.0);

for (int i=0; i<dim; i++){
  for (int j=0; j<dim; j++){
	  f_grad += constV(0.5*Kn1[i][j])*n1x[i]*n1x[j];
  }
}

// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double> > sfts1, sfts1c, sfts1cc;

for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  // Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c + b_p
	  sfts1[i][j] = constV(sfts_linear1[i][j])*c + constV(sfts_const1[i][j]);
	  sfts1c[i][j] = constV(sfts_linear1[i][j]);
	  sfts1cc[i][j] = constV(0.0);
  }
}

//compute E2=(E-E0)
dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V);

  }
}

//compute stress
//S=C*(E-E0)
dealii::VectorizedArray<double> CIJ_combined[2*dim-1+dim/3][2*dim-1+dim/3];

if (n_dependent_stiffness == true){
  for (unsigned int i=0; i<2*dim-1+dim/3; i++){
	  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
		  CIJ_combined[i][j] = this->userInputs.CIJ_list[0][i][j]*(constV(1.0)-h1V) + this->userInputs.CIJ_list[1][i][j]*h1V;
	  }
  }
  computeStress<dim>(CIJ_combined, E2, S);
}
else{
  computeStress<dim>(this->userInputs.CIJ_list[0], E2, S);
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




