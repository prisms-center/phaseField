// Definition of the variables in the model
#define num_var 3 // total number of variable 
#define variable_name {"c", "mu", "u"}
#define variable_type {"SCALAR","SCALAR","VECTOR"}
#define variable_eq_type {"PARABOLIC","PARABOLIC","ELLIPTIC"}
#define need_val {true, false, false}
#define need_grad {true, true, true}
#define need_hess {false, false, false} // Currently overridden based on value of "n_dependent_stiffness"
#define need_val_residual {true, true, false}
#define need_grad_residual {true, true, true} // ith term for ith Equ

// LHS denotes the CG-solver 
#define need_val_LHS {true, false, false}
#define need_grad_LHS {false, false, true}
#define need_hess_LHS {false, false, false}
#define need_val_residual_LHS {false, false, false}
#define need_grad_residual_LHS {false, false, true}

// Define Cahn-Hilliard parameters (no gradient energy terms)
#define McV 1.0 
#define KcV 0.5  

// Define Mechanical properties
#define c_dependent_stiffness true
#define n_dependent_stiffness true  // DELETE 
// Mechanical symmetry of the material and stiffness parameters
// Used throughout system if c_dependent_stiffness == false, used in c=0 phase if c_dependent_stiffness == true
#define MaterialModels {{"ANISOTROPIC"},{"ANISOTROPIC"}}
// #define MaterialConstantsV {1.0,1.0,0.385,0.5,0.0,0.0}
 #define MaterialConstants {{10.0,1.0,0.385,1.0,0.0,0.0},{10.0,1.0,0.385,1.0,0.0,0.0}}

//2D models:
//ISOTROPIC- (Plane Strain) 2 constants [E, nu]
//ANISOTROPIC- 6 constants [C11 C22 C33 C12 C13 C23]

// Used in n=1 phase if n_dependent_stiffness == true
//#define MaterialModelBetaV ANISOTROPIC
// #define MaterialConstantsBetaV {1.0,1.0,0.385,0.5,0.0,0.0}
//#define MaterialConstantsBetaV {10.0,1.0,0.385,1.0,0.0,0.0}

// Stress-free transformation strains
// Linear fits for the stress-free transformation strains in for sfts = sfts_linear * c + sfts_const
double sfts_const1[3][3] = {{0.00,0.00,0},{0,0.0,0},{0,0,0.00}};
double sfts_const2[3][3] = {{0.05,0.00,0},{0,0.01,0},{0,0,0.00}};

//define free energy expressions  
#define fV (0.25*c*c*(1.0-c)*(1.0-c)) // THIS IS USED JUST TO REPRESENT THE residuals below in a short form
#define fcV (0.5*(c*(1.0-c)*(1.0-c) -c*c*(1.0-c)))  // THIS IS USED JUST TO REPRESENT THE residuals below  in a short form

//define required residuals  
#define rmuV  (fcV -cDependentMisfitCH)
#define rmuxV (constV(KcV)*cx)
#define rcV   (c)
#define rcxV  (constV(-McV*timeStep)*mux)


template <int dim>
void generalizedProblem<dim>::residualRHS(const std::vector<modelVariable<dim>> & modelVariablesList,
												std::vector<modelResidual<dim>> & modelResidualsList,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//c
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

//mu
scalargradType mux = modelVariablesList[1].scalarGrad;

//u
vectorgradType ux = modelVariablesList[2].vectorGrad;
vectorgradType Rux; // COEFFICIENT OF GRADIENT OF w (TEST FUN FOR u)

// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double> > sfts, sftsc;

for (unsigned int i=0; i<dim; i++){  // LOOPING OVER THE 3X3 STRAIN MATRIX
for (unsigned int j=0; j<dim; j++){
	  // Polynomial fits for the stress-free transformation strains, of the form: sfts=eta0*I+(eta1-eta0)*I*h(c)=\epsilon0
	  // sfts = Stress Free Transformation Strain or Initial Strain
	  sfts[i][j]  = constV(sfts_const1[i][j]) + constV(sfts_const2[i][j])*c*c*(constV(3.0)-constV(2.0)*c);
	  sftsc[i][j] = constV(sfts_const2[i][j])*(constV(2.0)*c*(constV(3.0)-constV(2.0)*c) -constV(2.0)*c*c); //DERIVATIVE D\epsilon0_Dc
}
}

//compute E2=(E-E0)
dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim], SS[dim][dim];

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]) - sfts[i][j];
}
}

// compute stress
// S=C*(E-E0)
// Compute stress tensor (which is equal to the residual, Rux)
dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

if (c_dependent_stiffness == true){
for (unsigned int i=0; i<2*dim-1+dim/3; i++){
	  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
		  CIJ_combined[i][j] = CIJ_list[0][i][j]  + (c*c*(constV(3.0)-constV(2.0)*c))*(CIJ_list[1][i][j]-CIJ_list[0][i][j]); 
	  }  // CIJ_alpha => CIJ_0 Latex
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

// Compute DCijkl_Dc:(epsilon-epsilon0)
dealii::VectorizedArray<double> DCIJ_combined_Dc[CIJ_tensor_size][CIJ_tensor_size];

if (c_dependent_stiffness == true){
for (unsigned int i=0; i<2*dim-1+dim/3; i++){
	  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
		  DCIJ_combined_Dc[i][j] = (CIJ_list[1][i][j]-CIJ_list[0][i][j]) *(constV(2.0)*c*(constV(3.0)-constV(2.0)*c) -constV(2.0)*c*c); 
	  }  // CIJ_list[0] => CIJ_0 Latex
	}
computeStress<dim>(DCIJ_combined_Dc, E2, SS);
}

// Compute one of the stress terms in the chemical potential, cDependentMisfitCH = -C*(E-E0)*(D_E0_Dc) + (E-E0)*DC_Dc*(E-E0)
dealii::VectorizedArray<double> cDependentMisfitCH=constV(0.0);

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  cDependentMisfitCH+= - S[i][j]*(sftsc[i][j]) + SS[i][j] * constV(0.5)*E2[i][j];
}
}

modelResidualsList[0].scalarValueResidual = rcV;  // [*] this number * denotes the order of the dependent variables
modelResidualsList[0].scalarGradResidual = rcxV;

modelResidualsList[1].scalarValueResidual = rmuV;
modelResidualsList[1].scalarGradResidual = rmuxV;

modelResidualsList[2].vectorGradResidual = Rux;

}

// THIS IS FOR THE ELLIPTIC PART, I.E. THE MECH EQBM ... doesn't have the MISFIT PART OF STRAIN
template <int dim>
void generalizedProblem<dim>::residualLHS(const std::vector<modelVariable<dim>> & modelVarList,
		modelResidual<dim> & modelRes,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
//c
scalarvalueType c = modelVarList[0].scalarValue;

//u
vectorgradType ux = modelVarList[1].vectorGrad;
vectorgradType Rux;

// Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the dealii "symmetrize" function
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > E;
//E = symmetrize(ux); // Only works for Deal.II v8.3 and later
E = constV(0.5)*(ux + transpose(ux));

// Compute stress tensor (which is equal to the residual, Rux)
if (c_dependent_stiffness == true){
	dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > CIJ_combined;
		  CIJ_combined = CIJ_list[0]  + (c*c*(constV(3.0)-constV(2.0)*c))*(CIJ_list[1]-CIJ_list[0]); 
	computeStress<dim>(CIJ_combined, E, Rux);
}
else{
	computeStress<dim>(CIJ_list[0], E, Rux);
}

modelRes.vectorGradResidual = Rux;

}

template <int dim>
void generalizedProblem<dim>::energyDensity(const std::vector<modelVariable<dim>> & modelVarList,
											const dealii::VectorizedArray<double> & JxW_value,
											dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {
	scalarvalueType total_energy_density = constV(0.0);

//c
scalarvalueType c = modelVarList[0].scalarValue;
scalargradType cx = modelVarList[0].scalarGrad;

//mu
scalarvalueType mu = modelVarList[1].scalarValue;
scalargradType mux = modelVarList[1].scalarGrad;

//u
vectorgradType ux = modelVarList[2].vectorGrad;

scalarvalueType f_chem = fV;

scalarvalueType f_grad = constV(0.0); 

for (int i=0; i<dim; i++){
	  f_grad += constV(0.5)*constV(KcV)*cx[i]*cx[i];
}


// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double> > sfts, sftsc;

for (unsigned int i=0; i<dim; i++){  // LOOPING OVER THE 3X3 STRAIN MATRIX
for (unsigned int j=0; j<dim; j++){
	  // Polynomial fits for the stress-free transformation strains, of the form: sfts=eta0*I+(eta1-eta0)*I*h(c)=\epsilon0
	  // sfts = Stress Free Transformation Strain or Initial Strain
	  sfts[i][j]  = constV(sfts_const1[i][j]) + constV(sfts_const2[i][j])*c*c*(constV(3.0)-constV(2.0)*c);
}
}

//compute E2=(E-E0)
dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim]; 

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]) - sfts[i][j];
}
}

// compute stress
// S=C*(E-E0)
// Compute stress tensor (which is equal to the residual, Rux)
dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

if (c_dependent_stiffness == true){
for (unsigned int i=0; i<2*dim-1+dim/3; i++){
	  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
		  CIJ_combined[i][j] = CIJ_list[0][i][j]  + (c*c*(constV(3.0)-constV(2.0)*c))*(CIJ_list[1][i][j]-CIJ_list[0][i][j]); 
	  }  // CIJ_list[0] => CIJ_0 Latex
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




