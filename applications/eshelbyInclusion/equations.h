// List of variables and residual equations for the Eshelby inclusion example application

// =================================================================================
// Define the variables in the model
// =================================================================================
// The number of variables
#define num_var 1

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define variable_name {"u"}
#define variable_type {"VECTOR"}
#define variable_eq_type {"ELLIPTIC"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define need_val {false}
#define need_grad {true}
#define need_hess {false}

// Flags for whether the residual equation has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define need_val_residual {false}
#define need_grad_residual {true}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqn
// for the left-hand-side of the iterative solver for elliptic equations
#define need_val_LHS {false}
#define need_grad_LHS {true}
#define need_hess_LHS {false}

// Flags for whether the residual equation for the left-hand-side of the iterative
// solver for elliptic equations has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define need_val_residual_LHS {false}
#define need_grad_residual_LHS {true}

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".

// Define Mechanical properties
// Mechanical symmetry of the material and stiffness parameters
#define MaterialModels {"ISOTROPIC"}
#define MaterialConstants {{22.5,0.3}}

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
template <int dim>
void generalizedProblem<dim>::residualRHS(const std::vector<modelVariable<dim> > & modelVariablesList,
												std::vector<modelResidual<dim> > & modelResidualsList,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//u
vectorgradType ux = modelVariablesList[0].vectorGrad;
vectorgradType Rux;

dealii::VectorizedArray<double> sfts[dim][dim];

dealii::VectorizedArray<double> dist, a;

// Radius of the inclusion
a = constV(10.0);

// Distance from the center of the inclusion
dist = std::sqrt((q_point_loc[0]-constV(0.0))*(q_point_loc[0]-constV(0.0))
					+(q_point_loc[1]-constV(0.0))*(q_point_loc[1]-constV(0.0))
					+(q_point_loc[2]-constV(0.0))*(q_point_loc[2]-constV(0.0)));

// Calculation the stress-free transformation strain (the misfit)
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		if (i == j){

			sfts[i][j] = 0.01 * (0.5+ 0.5*( constV(1.0) - std::exp(-20.0*(dist-a)))/ (constV(1.0)+std::exp(-20.0*(dist-a))));

		}
		else {
			sfts[i][j] = 0.0;
		}
	}
}


//compute strain tensor
dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		E[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-sfts[i][j];
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
template <int dim>
void generalizedProblem<dim>::residualLHS(const std::vector<modelVariable<dim> > & modelVarList,
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
template <int dim>
void generalizedProblem<dim>::energyDensity(const std::vector<modelVariable<dim> > & modelVarList, const dealii::VectorizedArray<double> & JxW_value, dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {

	//u
	vectorgradType ux = modelVarList[0].vectorGrad;

	dealii::VectorizedArray<double> sfts[dim][dim];

	dealii::VectorizedArray<double> dist;

	dist = std::sqrt((q_point_loc[0]-constV(0.0))*(q_point_loc[0]-constV(0.0))
						+(q_point_loc[1]-constV(0.0))*(q_point_loc[1]-constV(0.0))
						+(q_point_loc[2]-constV(0.0))*(q_point_loc[2]-constV(0.0)));

	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			if (i == j){

				sfts[i][j] = 0.01 * (0.5+ 0.5*( constV(1.0) - std::exp(-20.0*(dist-constV(10.0))))/ (constV(1.0)+std::exp(-20.0*(dist-constV(10.0)))));

			}
			else {
				sfts[i][j] = 0.0;
			}
		}
	}


	//compute strain tensor
	dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			E[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-sfts[i][j];
		}
	}

	//compute stress tensor
	computeStress<dim>(CIJ_list[0], E, S);

	scalarvalueType f_el = constV(0.0);

	for (unsigned int i=0; i<dim; i++){
	  for (unsigned int j=0; j<dim; j++){
		  f_el += constV(0.5) * S[i][j]*E[i][j];
	  }
	}

	// Loop to step through each element of the vectorized arrays. Working with deal.ii
	// developers to see if there is a more elegant way to do this.
	assembler_lock.acquire ();
	for (unsigned i=0; i<f_el.n_array_elements;i++){
	  // For some reason, some of the values in this loop
	  if (f_el[i] > 1.0e-10){
		  this->energy+=f_el[i]*JxW_value[i];
	  }
	}
	assembler_lock.release ();


}




