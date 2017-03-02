// List of variables and residual equations for the diffusion example application

// =================================================================================
// Define the variables in the model
// =================================================================================
// The number of variables
#define num_var 1

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define variable_name {"c"}
#define variable_type {"SCALAR"}
#define variable_eq_type {"PARABOLIC"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define need_val {true}
#define need_grad {true}
#define need_hess  {false}

// Flags for whether the residual equation has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define need_val_residual {true}
#define need_grad_residual {true}

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".


// Diffusion constant
#define DcV 0.01

//define geometric expression for the source term
//defining a circle centered at (0.0,0.0) with radius spanX/10
#define exp1 (std::sqrt(x*x+y*y)<spanX/10)

//define temporal expression for the source term
#define exp2 (t<timeFinal/4)

// Residual equations
#define rcV   (c + timeStep*subdomain(exp1 && exp2, constV(1.0)) )
#define rcxV  (constV(-DcV*timeStep)*cx)

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


//c
scalarvalueType c = modelVariablesList[0].scalarValue;
scalargradType cx = modelVariablesList[0].scalarGrad;

double x=q_point_loc[0][0], y=q_point_loc[1][0];
double t=this->currentTime;

modelResidualsList[0].scalarValueResidual = rcV;
modelResidualsList[0].scalarGradResidual = rcxV;

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
void customPDE<dim,degree>::energyDensity(const std::vector<modelVariable<dim> > & modelVarList,
											const dealii::VectorizedArray<double> & JxW_value,
											dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {


}




