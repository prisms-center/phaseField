// Define required residuals

// ---------------------------------------------

template <int dim>
void generalizedProblem<dim>::residualRHS(const std::vector<modelVariable<dim>> & modelVariablesList,
												std::vector<modelResidual<dim>> & modelResidualsList,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//u
vectorgradType ux = modelVariablesList[0].vectorGrad;
vectorgradType Rux;

//compute strain tensor
dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		E[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]);
	}
}

//compute stress tensor
computeStress<dim>(CIJ, E, S);

//compute residual
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		Rux[i][j] = -S[i][j];
	}
}

modelResidualsList[0].vectorGradResidual = Rux;

}

template <int dim>
void generalizedProblem<dim>::residualLHS(const std::vector<modelVariable<dim>> & modelVarList,
		modelResidual<dim> & modelRes,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//u
vectorgradType ux = modelVarList[0].vectorGrad;
vectorgradType Rux;

//compute strain tensor
dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		E[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]);
	}
}

//compute stress tensor
computeStress<dim>(CIJ, E, S);

//compute residual
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		Rux[i][j] = S[i][j];
	}
}

modelRes.vectorGradResidual = Rux;

}

template <int dim>
void generalizedProblem<dim>::energyDensity(const std::vector<modelVariable<dim>> & modelVarList, const dealii::VectorizedArray<double> & JxW_value) {

}




