// Definition of the variables in the model
#define num_var 10
#define variable_name {"n1","n2","n3","n4","n5","n6","n7","n8","n9","n10"}
#define variable_type {"SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR"}
#define variable_eq_type {"PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC","PARABOLIC"}
#define need_val {true,true,true,true,true,true,true,true,true,true}
#define need_grad {true,true,true,true,true,true,true,true,true,true}
#define need_hess  {false,false,false,false,false,false,false,false,false,false}
#define need_val_residual {true,true,true,true,true,true,true,true,true,true}
#define need_grad_residual {true,true,true,true,true,true,true,true,true,true}

//define Allen-Cahn parameters
#define MnV 1.0
#define KnV 0.6
#define alpha 1.5

// All residuals defined below

// ---------------------------------------------

template <int dim>
void generalizedProblem<dim>::residualRHS(const std::vector<modelVariable<dim>> & modelVariablesList,
												std::vector<modelResidual<dim>> & modelResidualsList,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

dealii::VectorizedArray<double> fnV = constV(0.0);
scalargradType nx;

for (unsigned int i=0; i<num_var; i++){
	fnV = - modelVariablesList[i].scalarValue + modelVariablesList[i].scalarValue*modelVariablesList[i].scalarValue*modelVariablesList[i].scalarValue;
	nx = modelVariablesList[i].scalarGrad;
	for (unsigned int j=0; j<num_var; j++){
		if (i != j){
			fnV += constV(2.0*alpha) * modelVariablesList[i].scalarValue * modelVariablesList[j].scalarValue*modelVariablesList[j].scalarValue;
		}
	}
	modelResidualsList[i].scalarValueResidual = modelVariablesList[i].scalarValue-constV(timeStep*MnV)*fnV;
	modelResidualsList[i].scalarGradResidual = constV(-timeStep*KnV*MnV)*nx;
}


}

template <int dim>
void generalizedProblem<dim>::residualLHS(const std::vector<modelVariable<dim>> & modelVarList,
		modelResidual<dim> & modelRes,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}

template <int dim>
void generalizedProblem<dim>::energyDensity(const std::vector<modelVariable<dim>> & modelVarList,
											const dealii::VectorizedArray<double> & JxW_value,
											dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {
	scalarvalueType total_energy_density = constV(0.0);


}




