// Define required residuals
#define rcV   (c)
#define rcxTemp ( cx*((1.0-h1V-h2V-h3V)*faccV+(h1V+h2V+h3V)*fbccV) + n1x*((fbcV-facV)*hn1V) + n2x*((fbcV-facV)*hn2V) + n3x*((fbcV-facV)*hn3V) + grad_mu_el)
#define rcxV  (constV(-timeStep*McV)*rcxTemp)

#define rn1V   (n1-constV(timeStep*Mn1V)*((fbV-faV)*hn1V+nDependentMisfitAC1+heterMechAC1))
#define rn2V   (n2-constV(timeStep*Mn2V)*((fbV-faV)*hn2V+nDependentMisfitAC2+heterMechAC2))
#define rn3V   (n3-constV(timeStep*Mn3V)*((fbV-faV)*hn3V+nDependentMisfitAC3+heterMechAC3))
#define rn1xV  (constV(-timeStep*Mn1V)*Knx1)
#define rn2xV  (constV(-timeStep*Mn2V)*Knx2)
#define rn3xV  (constV(-timeStep*Mn3V)*Knx3)

// ---------------------------------------------

template <int dim>
void CoupledCHACMechanicsProblem<dim>::residualRHS(const std::vector<modelVariable<dim>> & modelVariablesList, std::vector<modelResidual<dim>> & modelResidualsList) const {

//c
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

//n1
scalarvalueType n1 = modelVariablesList[1].scalarValue;
scalargradType n1x = modelVariablesList[1].scalarGrad;

//n2
scalarvalueType n2 = modelVariablesList[2].scalarValue;
scalargradType n2x = modelVariablesList[2].scalarGrad;


//n3
scalarvalueType n3 = modelVariablesList[3].scalarValue;
scalargradType n3x = modelVariablesList[3].scalarGrad;

//u
vectorgradType ux = modelVariablesList[4].vectorGrad;
vectorgradType Rux;

vectorhessType uxx;

if (c_dependent_misfit == true){
	uxx = modelVariablesList[4].vectorHess;
}

// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double> > sfts1, sfts1c, sfts1cc, sfts2, sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc;

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  // Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c + b_p
	  sfts1[i][j] = constV(sfts_linear1[i][j])*c + constV(sfts_const1[i][j]);
	  sfts1c[i][j] = constV(sfts_linear1[i][j]);
	  sfts1cc[i][j] = constV(0.0);

	  // Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c + b_p
	  sfts2[i][j] = constV(sfts_linear2[i][j])*c + constV(sfts_const2[i][j]);
	  sfts2c[i][j] = constV(sfts_linear1[i][j]);
	  sfts2cc[i][j] = constV(0.0);

	  // Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c + b_p
	  sfts3[i][j] = constV(sfts_linear3[i][j])*c + constV(sfts_const3[i][j]);
	  sfts3c[i][j] = constV(sfts_linear3[i][j]);
	  sfts3cc[i][j] = constV(0.0);
}
}

//compute E2=(E-E0)
dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1strainV + sfts2[i][j]*h2strainV + sfts3[i][j]*h3strainV);

}
}

//compute stress
//S=C*(E-E0)
// Compute stress tensor (which is equal to the residual, Rux)
dealii::VectorizedArray<double> CIJ_combined[2*dim-1+dim/3][2*dim-1+dim/3];

if (n_dependent_stiffness == true){
dealii::VectorizedArray<double> sum_hV;
sum_hV = h1V+h2V+h3V;
for (unsigned int i=0; i<2*dim-1+dim/3; i++){
	  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
		  CIJ_combined[i][j] = constV(CIJ_alpha(i,j))*(constV(1.0)-sum_hV) + constV(CIJ_beta(i,j))*sum_hV;
	  }
}
computeStress<dim>(CIJ_combined, E2, S);
}
else{
computeStress<dim>(CIJ, E2, S);
}

// Fill residual corresponding to mechanics
// R=-C*(E-E0)

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  Rux[i][j] = - S[i][j];
}
}

// Compute one of the stress terms in the order parameter chemical potential, nDependentMisfitACp = C*(E-E0)*(E0_p*Hn)
dealii::VectorizedArray<double> nDependentMisfitAC1=constV(0.0);
dealii::VectorizedArray<double> nDependentMisfitAC2=constV(0.0);
dealii::VectorizedArray<double> nDependentMisfitAC3=constV(0.0);

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  nDependentMisfitAC1+=S[i][j]*(sfts1[i][j]);
	  nDependentMisfitAC2+=S[i][j]*(sfts2[i][j]);
	  nDependentMisfitAC3+=S[i][j]*(sfts3[i][j]);
}
}

nDependentMisfitAC1*=-hn1strainV;
nDependentMisfitAC2*=-hn2strainV;
nDependentMisfitAC3*=-hn3strainV;

// Compute the other stress term in the order parameter chemical potential, heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
dealii::VectorizedArray<double> heterMechAC1=constV(0.0);
dealii::VectorizedArray<double> heterMechAC2=constV(0.0);
dealii::VectorizedArray<double> heterMechAC3=constV(0.0);
dealii::VectorizedArray<double> S2[dim][dim];

if (n_dependent_stiffness == true){
computeStress<dim>(CIJ_diff, E2, S2);
for (unsigned int i=0; i<dim; i++){
	  for (unsigned int j=0; j<dim; j++){
		  heterMechAC1 += S2[i][j]*E2[i][j];
	  }
}
// Aside from HnpV, heterMechAC1, heterMechAC2, and heterMechAC3 are equal
heterMechAC2 = 0.5*hn2V*heterMechAC1;
heterMechAC3 = 0.5*hn3V*heterMechAC1;

heterMechAC1 = 0.5*hn1V*heterMechAC1;
}

// compute the stress term in the gradient of the concentration chemical potential, grad_mu_el = [C*(E-E0)*E0c]x, must be a vector with length dim
scalargradType grad_mu_el;

if (c_dependent_misfit == true){
	dealii::VectorizedArray<double> E3[dim][dim], S3[dim][dim];

	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			E3[i][j] =  -( sfts1c[i][j]*h1strainV + sfts2c[i][j]*h2strainV + sfts3c[i][j]*h3strainV);
		}
	}

	if (n_dependent_stiffness == true){
		computeStress<dim>(CIJ_combined, E3, S3);
	}
	else{
		computeStress<dim>(CIJ, E3, S3);
	}

	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			for (unsigned int k=0; k<dim; k++){
				grad_mu_el[k] += S3[i][j] * (constV(0.5)*(uxx[i][j][k]+uxx[j][i][k]) + E3[i][j]*cx[k]
														  - (sfts1[i][j]*hn1strainV*n1x[k] + sfts2[i][j]*hn2strainV*n2x[k] + sfts3[i][j]*hn3strainV*n3x[k]));

				grad_mu_el[k]+= - S[i][j] * (sfts1c[i][j]*hn1strainV*n1x[k] + sfts2c[i][j]*hn2strainV*n2x[k] + sfts3c[i][j]*hn3strainV*n3x[k]
														  + (sfts1cc[i][j]*h1strainV + sfts2cc[i][j]*h2strainV + sfts3cc[i][j]*h3strainV)*cx[k]);

				if (n_dependent_stiffness == true){
					grad_mu_el[k]+= - S2[i][j] * (sfts1c[i][j]*hn1V*n1x[k] + sfts2c[i][j]*hn2V*n2x[k] + sfts3c[i][j]*hn3V*n3x[k]);

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

modelResidualsList[1].scalarValueResidual = rn1V;
modelResidualsList[1].scalarGradResidual = rn1xV;

modelResidualsList[2].scalarValueResidual = rn2V;
modelResidualsList[2].scalarGradResidual = rn2xV;

modelResidualsList[3].scalarValueResidual = rn3V;
modelResidualsList[3].scalarGradResidual = rn3xV;

modelResidualsList[4].vectorGradResidual = Rux;

}

template <int dim>
void CoupledCHACMechanicsProblem<dim>::residualLHS(const std::vector<modelVariable<dim>> & modelVariablesList, std::vector<modelResidual<dim>> & modelResidualsList) const {

//n1
scalarvalueType n1 = modelVariablesList[0].scalarValue;

//n2
scalarvalueType n2 = modelVariablesList[1].scalarValue;


//n3
scalarvalueType n3 = modelVariablesList[2].scalarValue;

//u
vectorgradType ux = modelVariablesList[3].vectorGrad;
vectorgradType Rux;



// Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the dealii "symmetrize" function
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > E;
//E = symmetrize(ux); // Only works for Deal.II v8.3 and later
E = constV(0.5)*(ux + transpose(ux));

// Compute stress tensor (which is equal to the residual, Rux)
if (n_dependent_stiffness == true){
	dealii::Tensor<2, 2*dim-1+dim/3, dealii::VectorizedArray<double> > CIJ_combined;
	CIJ_combined = CIJ_alpha_tensor*(constV(1.0)-h1V-h2V-h3V) + CIJ_beta_tensor*(h1V+h2V+h3V);

	computeStress<dim>(CIJ_combined, E, Rux);
}
else{
	computeStress<dim>(CIJ, E, Rux);
}

modelResidualsList[0].vectorGradResidual = Rux;

}




