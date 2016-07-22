// Define Cahn-Hilliard parameters (no gradient energy terms)
#define McV 1.0

// Define Allen-Cahn parameters
#define MnV 150.0
#define KnV 0.5

//define free energy expressions
#define faV (-1.6704-4.776*c+5.1622*c*c-2.7375*c*c*c+1.3687*c*c*c*c)
#define facV (-4.776 + 10.3244*c - 8.2125*c*c + 5.4748*c*c*c)
#define faccV (10.3244-16.425*c+16.4244*c*c)
#define fbV (5.0*c*c-5.9746*c-1.5924)
#define fbcV (10.0*c-5.9746)
#define fbccV (10.0)
#define hV (10.0*n*n*n-15.0*n*n*n*n+6.0*n*n*n*n*n)
#define hnV (30.0*n*n-60.0*n*n*n+30.0*n*n*n*n)

//define required residuals
#define muxV ( cx*((1.0-hV)*faccV+hV*fbccV) + nx*((fbcV-facV)*hnV) )
#define rcV   (c)
#define rcxV  (constV(-McV*timeStep)*muxV)
#define rnV  (n-constV(timeStep*MnV)*(fbV-faV)*hnV)
#define rnxV (constV(-timeStep*KnV*MnV)*nx)


// ---------------------------------------------

template <int dim>
void generalizedProblem<dim>::residualRHS(const std::vector<modelVariable<dim>> & modelVariablesList, std::vector<modelResidual<dim>> & modelResidualsList) const {

//c
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

//n
scalarvalueType n = modelVariablesList[1].scalarValue;
scalargradType nx = modelVariablesList[1].scalarGrad;

modelResidualsList[0].scalarValueResidual = rcV;
modelResidualsList[0].scalarGradResidual = rcxV;

modelResidualsList[1].scalarValueResidual = rnV;
modelResidualsList[1].scalarGradResidual = rnxV;

}

template <int dim>
void generalizedProblem<dim>::residualLHS(const std::vector<modelVariable<dim>> & modelVarList, modelResidual<dim> & modelRes) const {

}

template <int dim>
void generalizedProblem<dim>::energyDensity(const std::vector<modelVariable<dim>> & modelVarList, const dealii::VectorizedArray<double> & JxW_value) {
	scalarvalueType total_energy_density = constV(0.0);

//c
scalarvalueType c = modelVarList[0].scalarValue;
scalargradType cx = modelVarList[0].scalarGrad;

//n1
scalarvalueType n = modelVarList[1].scalarValue;
scalargradType nx = modelVarList[1].scalarGrad;

scalarvalueType f_chem = (constV(1.0)-hV)*faV + hV*fbV;

scalarvalueType f_grad = constV(0.5*KnV)*nx*nx;

total_energy_density = f_chem + f_grad;

assembler_lock.acquire ();
for (unsigned i=0; i<c.n_array_elements;i++){
  // For some reason, some of the values in this loop
  if (c[i] > 1.0e-10){
	  this->energy+=total_energy_density[i]*JxW_value[i];
	  this->energy_components[0]+= f_chem[i]*JxW_value[i];
	  this->energy_components[1]+= f_grad[i]*JxW_value[i];
  }
}
assembler_lock.release ();
}




