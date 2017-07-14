// List of residual equations for the coupled Allen-Cahn example application

/*
// =================================================================================
// setVariableAttributes
// =================================================================================
// [INSERT EXPLANATION HERE]
template <int dim, int degree>
void customPDE<dim,degree>::setVariableAttributes(){

// Variable 0
userInputs.set_variable_name(0,”n”);
userInputs.set_variable_type(0,SCALAR);
userInputs.set_variable_equation_type(0,PARABOLIC);

userInputs.need_variable_value(0,true);
userInputs.need_variable_gradient(0,true);
userInputs.need_variable_hessian(0,false);

userInputs.need_value_residual_term(0,true);
userInputs.need_gradient_residual_term(0,true);

}
*/

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
void customPDE<dim,degree>::residualRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// The order parameter and its derivatives (names here should match those in the macros above)
scalarvalueType n = variable_list.get_scalar_value(0); //modelVariablesList[0].scalarValue;
scalargradType nx = variable_list.get_scalar_gradient(0); //modelVariablesList[0].scalarGrad;

// Parameters in the residual equations and expressions for the residual equations
// can be set here.
scalarvalueType fnV = (4.0*n*(n-1.0)*(n-0.5));
scalarvalueType rnV = (n-constV(userInputs.dtValue*MnV)*fnV);
scalargradType rnxV = (-constV(userInputs.dtValue*KnV*MnV)*nx);

// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
variable_list.set_scalar_value_residual_term(0,rnV);
variable_list.set_scalar_gradient_residual_term(0,rnxV);

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
void customPDE<dim,degree>::residualLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}
