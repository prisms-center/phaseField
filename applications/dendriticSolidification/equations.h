// List of variables and residual equations for the coupled Allen-Cahn/Cahn-Hilliard example application

// =================================================================================
// Define the variables in the model
// =================================================================================
// The number of variables
#define num_var 3

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define variable_name {"T", "n", "mu"}
#define variable_type {"SCALAR","SCALAR","SCALAR"}
#define variable_eq_type {"PARABOLIC","PARABOLIC","PARABOLIC"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define need_val {true, true, true}
#define need_grad {true, true, false}
#define need_hess {false, false, false}

// Flags for whether the residual equation has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define need_val_residual {true, true, true}
#define need_grad_residual {true, false, true}

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".

// Heat diffusion constant
#define DV 10.0

// Free energy expression (or rather, its derivative)
#define lambdaV (DV/0.6267)
#define fnV (-(n-constV(lambdaV)*T*(constV(1.0)-n*n))*(constV(1.0)-n*n))

// Undercooling
#define deltaV (0.05*8.0)

// Anisotropy paramter
#define epsilonM (0.025*5.0)

// Rotation angle
#define theta0 0.0

// Symmetry factor
#define mult 4.0

// Isotropic gradient energy coefficient
#define W0 0.1

// Anisotropic gradient energy coefficient, its derivative and square
#define WV (constV(W0)*(constV(1.0)+constV(epsilonM)*std::cos(constV(mult)*theta-constV(theta0))))
#define WTV (constV(W0)*(-constV(epsilonM)*constV(mult)*std::sin(constV(mult)*theta-constV(theta0))))
#define tauV (WV*WV)

// Define required residuals (theta and aniso defined in residualRHS)
#define rTV   (T-constV(0.5)*mu*constV(timeStep)/tauV)
#define rTxV  (constV(-DV*timeStep)*Tx)
#define rnV  (n-constV(timeStep)*mu/tauV)
#define rmuV (fnV)
#define rmuxV (-aniso)



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
void customPDE<dim,degree>::residualRHS(const std::vector<modelVariable<dim> > & modelVariablesList,
												std::vector<modelResidual<dim> > & modelResidualsList,
												dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// The temperature and its derivatives (names here should match those in the macros above)
scalarvalueType T = modelVariablesList[0].scalarValue;
scalargradType Tx = modelVariablesList[0].scalarGrad;

// The order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n = modelVariablesList[1].scalarValue;
scalargradType nx = modelVariablesList[1].scalarGrad;

// The order parameter chemical potential and its derivatives (names here should match those in the macros above)
scalarvalueType mu = modelVariablesList[2].scalarValue;


scalarvalueType theta;

for (unsigned i=0; i< n.n_array_elements;i++){
	theta[i] = std::atan2(nx[1][i],nx[0][i]);
}

scalargradType aniso;
aniso[0] = -WV*WV*nx[0]+WV*WTV*nx[1];
aniso[1] = -WV*WV*nx[1]-WV*WTV*nx[0];


// Residuals for the equation to evolve the concentration (names here should match those in the macros above)
modelResidualsList[0].scalarValueResidual = rTV;
modelResidualsList[0].scalarGradResidual = rTxV;

// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
modelResidualsList[1].scalarValueResidual = rnV;

// Residuals for the equation to evolve the order parameter chemical potential (names here should match those in the macros above)
modelResidualsList[2].scalarValueResidual = rmuV;
modelResidualsList[2].scalarGradResidual = rmuxV;

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
void customPDE<dim,degree>::residualLHS(const std::vector<modelVariable<dim> > & modelVarList,
		modelResidual<dim> & modelRes,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
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
void customPDE<dim,degree>::energyDensity(const std::vector<modelVariable<dim> > & modelVariablesList,
											const dealii::VectorizedArray<double> & JxW_value,
											dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {


}




