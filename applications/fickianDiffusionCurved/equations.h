// List of variables and residual equations for the diffusion example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"c");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,PARABOLIC);

	set_need_value					(0,true);
	set_need_gradient				(0,true);
	set_need_hessian				(0,false);

	set_need_value_residual_term	(0,true);
	set_need_gradient_residual_term	(0,true);

}

// =================================================================================
// residualRHS
// =================================================================================
// This function calculates the residual equations for each variable. It takes
// "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of
// that quadrature point is given by "q_point_loc". The function outputs residuals
// to variable_list. The index for each variable in this list corresponds to
// the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::residualRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {


//c
scalarvalueType c = variable_list.get_scalar_value(0);
scalargradType cx = variable_list.get_scalar_gradient(0);

scalarvalueType x=q_point_loc[0], y=q_point_loc[1];
double t=this->currentTime;
double T = this->userInputs.finalTime;

double t_1 = 0.2*T;
double tau_1 = 0.2*T;
scalarvalueType x_1 = constV(0.6*userInputs.domain_size[0]);
scalarvalueType y_1 = constV(0.2*userInputs.domain_size[1]);
scalarvalueType L_1 = constV(0.01*(userInputs.domain_size[0]+userInputs.domain_size[1]));

double t_2 = 0.6*T;
double tau_2 = 0.2*T;
scalarvalueType x_2 = constV(0.3*userInputs.domain_size[0]);
scalarvalueType y_2 = constV(0.7*userInputs.domain_size[1]);
scalarvalueType L_2 = constV(0.01*(userInputs.domain_size[0]+userInputs.domain_size[1]));

scalarvalueType source_term1 = 100.0*std::exp( - (t-t_1)/tau_1 * (t-t_1)/tau_1 )
								*std::exp( -((x-x_1)*(x-x_1)+(y-y_1)*(y-y_1))/(L_1*L_1) );

scalarvalueType source_term2 = 100.0*std::exp( - (t-t_2)/tau_2 * (t-t_2)/tau_2 )
								*std::exp( -((x-x_2)*(x-x_2)+(y-y_2)*(y-y_2))/(L_2*L_2) );


// Residual equations
scalarvalueType rcV = (c + userInputs.dtValue*(source_term1 + source_term2) );
scalargradType rcxV = (constV(-DcV*userInputs.dtValue)*cx);


variable_list.set_scalar_value_residual_term(0,rcV);
variable_list.set_scalar_gradient_residual_term(0,rcxV);

}

// =================================================================================
// residualLHS (needed only if at least one equation is elliptic)
// =================================================================================
// This function calculates the residual equations for the iterative solver for
// elliptic equations.for each variable. It takes "variable_list" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is given
// by "q_point_loc". The function outputs residual terms to "variable_list"
// for the left-hand-side of the residual equation for the iterative solver. The
// index for each variable in this list corresponds to
// the index given at the top of this file. If there are multiple elliptic equations,
// conditional statements should be used to ensure that the correct residual is
// being submitted. The index of the field being solved can be accessed by
// "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::residualLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}
