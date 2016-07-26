// Definition of the variables in the model
#define num_var 5
#define variable_name {"c", "n1", "n2", "n3", "u"}
#define variable_type {"SCALAR","SCALAR","SCALAR","SCALAR","VECTOR"}
#define variable_eq_type {"PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","ELLIPTIC"}
#define need_val {true, true, true, true, false}
#define need_grad {true, true, true, true, true}
#define need_hess {false, false, false, false, false} // Currently overridden based on value of "n_dependent_stiffness"
#define need_val_residual {true, true, true, true, false}
#define need_grad_residual {true, true, true, true, true}

#define need_val_LHS {false, true, true, true, false}
#define need_grad_LHS {false, false, false, false, true}
#define need_hess_LHS {false, false, false, false, false}
#define need_val_residual_LHS {false, false, false, false, false}
#define need_grad_residual_LHS {false, false, false, false, true}

// Define Cahn-Hilliard parameters (no gradient energy terms)
#define McV 1.0

// Define Allen-Cahn parameters
#define Mn1V 100.0
#define Mn2V 100.0
#define Mn3V 100.0

double Kn1[3][3]={{0.03,0,0},{0,0.007,0},{0,0,1.0}};
double Kn2[3][3]={{0.01275,-0.009959,0},{-0.009959,0.02425,0},{0,0,1.0}};
double Kn3[3][3]={{0.01275,0.009959,0},{0.009959,0.02425,0},{0,0,1.0}};

// Define Mechanical properties
#define n_dependent_stiffness true
// Mechanical symmetry of the material and stiffness parameters
// Used throughout system if n_dependent_stiffness == false, used in n=0 phase if n_dependent_stiffness == true
//#define MaterialModelV ISOTROPIC
//#define MaterialConstantsV {2.0,0.3}

// Used in n=1 phase if n_dependent_stiffness == true
//#define MaterialModelBetaV ISOTROPIC
//#define MaterialConstantsBetaV {2.5,0.3}

#define MaterialModels {{"ISOTROPIC"},{"ISOTROPIC"}}
#define MaterialConstants {{2.0,0.3},{2.5,0.3}}

// Stress-free transformation strains
// Linear fits for the stress-free transformation strains in for sfts = sfts_linear * c + sfts_const
double sfts_linear1[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double sfts_const1[3][3] = {{0.0345,0,0},{0,0.0185,0},{0,0,-0.00270}};

double sfts_linear2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double sfts_const2[3][3]={{0.0225,-0.0069,0},{-0.0069,0.0305,0},{0,0,-0.00270}};

double sfts_linear3[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
double sfts_const3[3][3]={{0.0225, 0.0069,0},{0.0069,0.0305,0},{0,0,-0.00270}};

//define free energy expressions
#define faV (-1.6704-4.776*c+5.1622*c*c-2.7375*c*c*c+1.3687*c*c*c*c)
#define facV (-4.776 + 10.3244*c - 8.2125*c*c + 5.4748*c*c*c)
#define faccV (10.3244-16.425*c+16.4244*c*c)
#define fbV (5.0*c*c-5.9746*c-1.5924)
#define fbcV (10.0*c-5.9746)
#define fbccV (10.0)
#define h1V (10.0*n1*n1*n1-15.0*n1*n1*n1*n1+6.0*n1*n1*n1*n1*n1)
#define h2V (10.0*n2*n2*n2-15.0*n2*n2*n2*n2+6.0*n2*n2*n2*n2*n2)
#define h3V (10.0*n3*n3*n3-15.0*n3*n3*n3*n3+6.0*n3*n3*n3*n3*n3)
#define hn1V (30.0*n1*n1-60.0*n1*n1*n1+30.0*n1*n1*n1*n1)
#define hn2V (30.0*n2*n2-60.0*n2*n2*n2+30.0*n2*n2*n2*n2)
#define hn3V (30.0*n3*n3-60.0*n3*n3*n3+30.0*n3*n3*n3*n3)

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
void generalizedProblem<dim>::residualRHS(const std::vector<modelVariable<dim>> & modelVariablesList,
												std::vector<modelResidual<dim>> & modelResidualsList,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
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
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V + sfts2[i][j]*h2V + sfts3[i][j]*h3V);

}
}

//compute stress
//S=C*(E-E0)
// Compute stress tensor (which is equal to the residual, Rux)
dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

if (n_dependent_stiffness == true){
dealii::VectorizedArray<double> sum_hV;
sum_hV = h1V+h2V+h3V;
for (unsigned int i=0; i<2*dim-1+dim/3; i++){
	  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
		  CIJ_combined[i][j] = CIJ_list[0][i][j]*(constV(1.0)-sum_hV) + CIJ_list[1][i][j]*sum_hV;
	  }
}
computeStress<dim>(CIJ_combined, E2, S);
}
else{
computeStress<dim>(CIJ_list[0], E2, S);
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

nDependentMisfitAC1*=-hn1V;
nDependentMisfitAC2*=-hn2V;
nDependentMisfitAC3*=-hn3V;

// Compute the other stress term in the order parameter chemical potential, heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
dealii::VectorizedArray<double> heterMechAC1=constV(0.0);
dealii::VectorizedArray<double> heterMechAC2=constV(0.0);
dealii::VectorizedArray<double> heterMechAC3=constV(0.0);
dealii::VectorizedArray<double> S2[dim][dim];

if (n_dependent_stiffness == true){
//computeStress<dim>(CIJ_diff, E2, S2);
	computeStress<dim>(CIJ_list[1]-CIJ_list[0], E2, S2);
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
			E3[i][j] =  -( sfts1c[i][j]*h1V + sfts2c[i][j]*h2V + sfts3c[i][j]*h3V);
		}
	}

	if (n_dependent_stiffness == true){
		computeStress<dim>(CIJ_combined, E3, S3);
	}
	else{
		computeStress<dim>(CIJ_list[0], E3, S3);
	}

	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			for (unsigned int k=0; k<dim; k++){
				grad_mu_el[k] += S3[i][j] * (constV(0.5)*(uxx[i][j][k]+uxx[j][i][k]) + E3[i][j]*cx[k]
														  - (sfts1[i][j]*hn1V*n1x[k] + sfts2[i][j]*hn2V*n2x[k] + sfts3[i][j]*hn3V*n3x[k]));

				grad_mu_el[k]+= - S[i][j] * (sfts1c[i][j]*hn1V*n1x[k] + sfts2c[i][j]*hn2V*n2x[k] + sfts3c[i][j]*hn3V*n3x[k]
														  + (sfts1cc[i][j]*h1V + sfts2cc[i][j]*h2V + sfts3cc[i][j]*h3V)*cx[k]);

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
void generalizedProblem<dim>::residualLHS(const std::vector<modelVariable<dim>> & modelVarList,
		modelResidual<dim> & modelRes,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//n1
scalarvalueType n1 = modelVarList[0].scalarValue;

//n2
scalarvalueType n2 = modelVarList[1].scalarValue;


//n3
scalarvalueType n3 = modelVarList[2].scalarValue;

//u
vectorgradType ux = modelVarList[3].vectorGrad;
vectorgradType Rux;

// Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the dealii "symmetrize" function
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > E;
//E = symmetrize(ux); // Only works for Deal.II v8.3 and later
E = constV(0.5)*(ux + transpose(ux));

// Compute stress tensor (which is equal to the residual, Rux)
if (n_dependent_stiffness == true){
	dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > CIJ_combined;
	CIJ_combined = CIJ_list[0]*(constV(1.0)-h1V-h2V-h3V) + CIJ_list[1]*(h1V+h2V+h3V);

	computeStress<dim>(CIJ_combined, E, Rux);
}
else{
	computeStress<dim>(CIJ_list[0], E, Rux);
}

modelRes.vectorGradResidual = Rux;

}

template <int dim>
void generalizedProblem<dim>::energyDensity(const std::vector<modelVariable<dim>> & modelVarList, const dealii::VectorizedArray<double> & JxW_value) {
	scalarvalueType total_energy_density = constV(0.0);

//c
scalarvalueType c = modelVarList[0].scalarValue;
scalargradType cx = modelVarList[0].scalarGrad;

//n1
scalarvalueType n1 = modelVarList[1].scalarValue;
scalargradType n1x = modelVarList[1].scalarGrad;

//n2
scalarvalueType n2 = modelVarList[2].scalarValue;
scalargradType n2x = modelVarList[2].scalarGrad;


//n3
scalarvalueType n3 = modelVarList[3].scalarValue;
scalargradType n3x = modelVarList[3].scalarGrad;

//u
vectorgradType ux = modelVarList[4].vectorGrad;

scalarvalueType f_chem = (constV(1.0)-(h1V+h2V+h3V))*faV + (h1V+h2V+h3V)*fbV;

scalarvalueType f_grad = constV(0.0);

for (int i=0; i<dim; i++){
  for (int j=0; j<dim; j++){
	  f_grad += constV(0.5*Kn1[i][j])*n1x[i]*n1x[j];
  }
}
#if num_sop>1
for (int i=0; i<dim; i++){
  for (int j=0; j<dim; j++){
	  f_grad += constV(0.5*Kn2[i][j])*n2x[i]*n2x[j];
  }
}
#endif
#if num_sop>2
for (int i=0; i<dim; i++){
  for (int j=0; j<dim; j++){
	  f_grad += constV(0.5*Kn3[i][j])*n3x[i]*n3x[j];
  }
}
#endif


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
	  //E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V + sfts2[i][j]*h2V + sfts3[i][j]*h3V);
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V + sfts2[i][j]*h2V + sfts3[i][j]*h3V);

  }
}

//compute stress
//S=C*(E-E0)
dealii::VectorizedArray<double> CIJ_combined[2*dim-1+dim/3][2*dim-1+dim/3];

if (n_dependent_stiffness == true){
  dealii::VectorizedArray<double> sum_hV;
  sum_hV = h1V+h2V+h3V;
  for (unsigned int i=0; i<2*dim-1+dim/3; i++){
	  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
		  CIJ_combined[i][j] = CIJ_list[0][i][j]*(constV(1.0)-sum_hV) + CIJ_list[1][i][j]*sum_hV;
	  }
  }
  computeStress<dim>(CIJ_combined, E2, S);
}
else{
  computeStress<dim>(CIJ_list[0], E2, S);
}

scalarvalueType f_el = constV(0.0);

for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  f_el += constV(0.5) * S[i][j]*E2[i][j];
  }
}

total_energy_density = f_chem + f_grad + f_el;

assembler_lock.acquire ();
for (unsigned i=0; i<c.n_array_elements;i++){
  // For some reason, some of the values in this loop
  if (c[i] > 1.0e-10){
	  this->energy+=total_energy_density[i]*JxW_value[i];
	  this->energy_components[0]+= f_chem[i]*JxW_value[i];
	  this->energy_components[1]+= f_grad[i]*JxW_value[i];
	  this->energy_components[2]+= f_el[i]*JxW_value[i];
  }
}
assembler_lock.release ();
}




