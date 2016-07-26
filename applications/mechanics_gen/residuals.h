// Define required residuals
// Definition of the variables in the model
#define num_var 1
#define variable_name {"u"}
#define variable_type {"VECTOR"}
#define variable_eq_type {"ELLIPTIC"}
#define need_val {false}
#define need_grad {true}
#define need_hess {false}
#define need_val_residual {false}
#define need_grad_residual {true}

#define need_val_LHS {false}
#define need_grad_LHS {true}
#define need_hess_LHS {false}
#define need_val_residual_LHS {false}
#define need_grad_residual_LHS {true}


// Define Mechanical properties
// Mechanical symmetry of the material and stiffness parameters
#define MaterialModels {"ISOTROPIC"}
#define MaterialConstants {{2.0,0.3}}


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
computeStress<dim>(CIJ_list[0], E, S);

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
computeStress<dim>(CIJ_list[0], E, S);

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




