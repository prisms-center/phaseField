// Definition of the variables in the model
#define num_var 2
#define variable_name {"c", "mu"}
#define variable_type {"SCALAR", "SCALAR"}
#define variable_eq_type {"PARABOLIC", "PARABOLIC"}
#define need_val {true, true}
#define need_grad {true, true}
#define need_hess  {false, false}
#define need_val_residual {true, true}
#define need_grad_residual {true, true}

//define Cahn-Hilliard parameters
#define McV 1.0
#define KcV 1.0

//define the free energy and its derivative with respect to n
#define fV (c*c*c*c - 2.0*c*c*c + c*c)
#define fcV (4.0*c*(c-1.0)*(c-0.5))

//define required residuals
#define rmuV  (fcV)
#define rmuxV (constV(KcV)*cx)
#define rcV   (c)
#define rcxV  (constV(-McV*timeStep)*mux)


// ---------------------------------------------

template <int dim>
void generalizedProblem<dim>::residualRHS(const std::vector<modelVariable<dim>> & modelVariablesList,
												std::vector<modelResidual<dim>> & modelResidualsList,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

	// c
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

// mu
scalarvalueType mu = modelVariablesList[1].scalarValue;
scalargradType mux = modelVariablesList[1].scalarGrad;



modelResidualsList[0].scalarValueResidual = rcV;
modelResidualsList[0].scalarGradResidual = rcxV;

modelResidualsList[1].scalarValueResidual = rmuV;
modelResidualsList[1].scalarGradResidual = rmuxV;

}

template <int dim>
void generalizedProblem<dim>::residualLHS(const std::vector<modelVariable<dim>> & modelVarList,
		modelResidual<dim> & modelRes,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}

template <int dim>
void generalizedProblem<dim>::energyDensity(const std::vector<modelVariable<dim>> & modelVariablesList, const dealii::VectorizedArray<double> & JxW_value) {
	scalarvalueType total_energy_density = constV(0.0);

// c
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

// mu
scalarvalueType mu = modelVariablesList[1].scalarValue;
scalargradType mux = modelVariablesList[1].scalarGrad;


scalarvalueType f_chem = fV;

scalarvalueType f_grad = constV(0.0);

for (int i=0; i<dim; i++){
  for (int j=0; j<dim; j++){
	  f_grad += constV(0.5*KcV)*cx[i]*cx[j];
  }
}

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




