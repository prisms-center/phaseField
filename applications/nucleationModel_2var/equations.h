
// List of variables and residual equations for the coupled Allen-Cahn/Cahn-Hilliard example application

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

	set_allowed_to_nucleate			(0, false);
	set_need_value_nucleation		(0, true);

    // Variable 1
	set_variable_name				(1,"n1");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,PARABOLIC);

	set_need_value					(1,true);
	set_need_gradient				(1,true);
	set_need_hessian				(1,false);

	set_need_value_residual_term	(1,true);
	set_need_gradient_residual_term	(1,true);

	set_allowed_to_nucleate			(1, true);
	set_need_value_nucleation		(1, true);

	// Variable 1
	set_variable_name				(2,"n2");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,PARABOLIC);

	set_need_value					(2,true);
	set_need_gradient				(2,true);
	set_need_hessian				(2,false);

	set_need_value_residual_term	(2,true);
	set_need_gradient_residual_term	(2,true);

	set_allowed_to_nucleate			(2, true);
	set_need_value_nucleation		(2, true);
}

// =================================================================================
// Define the Model and residual equations
// =================================================================================
// The model free energy equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".

// Free energy for each phase and their first and second derivatives
#define faV (A0+A2*(c_alpha-calmin)*(c_alpha-calmin))
#define facV (2.0*A2*(c_alpha-calmin))
#define faccV (2.0*A2)
#define fbV (B0+B2*(c_beta-cbtmin)*(c_beta-cbtmin))
#define fbcV (2.0*B2*(c_beta-cbtmin))
#define fbccV (2.0*B2)

// Interpolation function and its derivative
#define h1V (3.0*n1*n1 - 2.0*n1*n1*n1)
#define h1nV (6.0*n1 - 6.0*n1*n1)

#define h2V (3.0*n2*n2 - 2.0*n2*n2*n2)
#define h2nV (6.0*n2 - 6.0*n2*n2)

// KKS model c_alpha and c_beta as a function of c and h
#define c_alpha ((B2*(c-cbtmin*(h1V+h2V)) + A2*calmin*(h1V+h2V))/(A2*(h1V+h2V)+B2*(1.0-(h1V+h2V))))
#define c_beta ((A2*(c-calmin*(1.0-(h1V+h2V)))+B2*cbtmin*(1.0-(h1V+h2V)))/(A2*(h1V+h2V)+B2*(1.0-(h1V+h2V))))

// Double-Well function (can be used to tune the interfacial energy)
#define fbarrierV (n1*n1 - 2.0*n1*n1*n1 + n1*n1*n1*n1 + n2*n2 - 2.0*n2*n2*n2 + n2*n2*n2*n2)
#define fbarriern1V (2.0*n1 - 6.0*n1*n1 + 4.0*n1*n1*n1)
#define fbarriern2V (2.0*n2 - 6.0*n2*n2 + 4.0*n2*n2*n2)

// Residual equations
// For concentration
#define term_muxV (cx + (c_alpha - c_beta)*(h1nV*n1x + h2nV*n2x))
#define rcV   (c)
#define rcxV  (constV(-McV*userInputs.dtValue)*term_muxV)
//For order parameter (gamma is a variable order parameter mobility factor)
#define rn1V   (n1-constV(userInputs.dtValue*MnV)*gamma*((fbV-faV)*h1nV - (c_beta-c_alpha)*fbcV*h1nV + W_barrier*fbarriern1V))
#define rn1xV  (constV(-userInputs.dtValue*KnV*MnV)*gamma*n1x)
#define rn2V   (n2-constV(userInputs.dtValue*MnV)*gamma*((fbV-faV)*h2nV - (c_beta-c_alpha)*fbcV*h2nV + W_barrier*fbarriern2V))
#define rn2xV  (constV(-userInputs.dtValue*KnV*MnV)*gamma*n2x)

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

// The concentration and its derivatives
scalarvalueType c = variable_list.get_scalar_value(0);
scalargradType cx = variable_list.get_scalar_gradient(0);

// The order parameter and its derivatives
scalarvalueType n1 = variable_list.get_scalar_value(1);
scalargradType n1x = variable_list.get_scalar_gradient(1);

// The order parameter and its derivatives
scalarvalueType n2 = variable_list.get_scalar_value(2);
scalargradType n2x = variable_list.get_scalar_gradient(2);

// -------------------------------------------------
// Nucleation expressions
// -------------------------------------------------
dealii::VectorizedArray<double> source_term1 = constV(0.0);
dealii::VectorizedArray<double> source_term2 = constV(0.0);
dealii::VectorizedArray<double> gamma = constV(1.0);
seedNucleus(q_point_loc,source_term1,source_term2,gamma);
// -------------------------------------------------

// Residuals for the equation to evolve the concentration
variable_list.set_scalar_value_residual_term(0,rcV);
variable_list.set_scalar_gradient_residual_term(0,rcxV);

// Residuals for the equation to evolve the order parameter
variable_list.set_scalar_value_residual_term(1,rn1V+source_term1);
variable_list.set_scalar_gradient_residual_term(1,rn1xV);

// Residuals for the equation to evolve the order parameter
variable_list.set_scalar_value_residual_term(2,rn2V+source_term2);
variable_list.set_scalar_gradient_residual_term(2,rn2xV);

}

// =================================================================================
// seedNucleus: a function particular to this app
// =================================================================================
template <int dim,int degree>
void customPDE<dim,degree>::seedNucleus(const dealii::Point<dim, dealii::VectorizedArray<double> > & q_point_loc,
	dealii::VectorizedArray<double> & source_term1,
	dealii::VectorizedArray<double> & source_term2,
	dealii::VectorizedArray<double> & gamma) const {

	for (typename std::vector<nucleus<dim> >::const_iterator thisNucleus=this->nuclei.begin(); thisNucleus!=this->nuclei.end(); ++thisNucleus){

		if (thisNucleus->seededTime + thisNucleus->seedingTime > this->currentTime){

			// Calculate the weighted distance function to the order parameter freeze boundary (weighted_dist = 1.0 on that boundary)
			dealii::VectorizedArray<double> weighted_dist = constV(0.0);
			dealii::Tensor<1,dim,dealii::VectorizedArray<double> > shortest_edist_tensor;
			for (unsigned int j=0; j<dim; j++){
				shortest_edist_tensor[j] = thisNucleus->center(j) - q_point_loc(j);

				if (userInputs.BC_list[userInputs.nucleation_parameters_list_index.at(thisNucleus->orderParameterIndex)].var_BC_type[2*j]==PERIODIC){
					for (unsigned k=0; k<gamma.n_array_elements;k++){
						shortest_edist_tensor[j][k] = shortest_edist_tensor[j][k]-round(shortest_edist_tensor[j][k]/userInputs.domain_size[j])*userInputs.domain_size[j];
					}
				}
			}
			shortest_edist_tensor = userInputs.nucleation_parameters_list[userInputs.nucleation_parameters_list_index.at(thisNucleus->orderParameterIndex)].rotation_matrix * shortest_edist_tensor;
			for (unsigned int j=0; j<dim; j++){
				shortest_edist_tensor[j] /= constV(userInputs.nucleation_parameters_list[userInputs.nucleation_parameters_list_index.at(thisNucleus->orderParameterIndex)].freeze_semiaxes[j]);
			}
			weighted_dist = shortest_edist_tensor.norm_square();


			for (unsigned i=0; i<gamma.n_array_elements;i++){
				if (weighted_dist[i] <= 1.0){
					gamma[i] = 0.0;

					// Seed a nucleus if it was added to the list of nuclei this time step
					if (thisNucleus->seedingTimestep == this->currentIncrement){
						// Find the weighted distance to the outer edge of the nucleus and use it to calculate the order parameter source term
						double r = 0.0;
						double avg_semiaxis = 0.0;
						dealii::Tensor<1,dim,double> shortest_edist_tensor;
						for (unsigned int j=0; j<dim; j++){
							shortest_edist_tensor[j] = thisNucleus->center(j) - q_point_loc(j)[i];

							if (userInputs.BC_list[userInputs.nucleation_parameters_list_index.at(thisNucleus->orderParameterIndex)].var_BC_type[2*j]==PERIODIC){
								shortest_edist_tensor[j] = shortest_edist_tensor[j]-round(shortest_edist_tensor[j]/userInputs.domain_size[j])*userInputs.domain_size[j];
							}
							avg_semiaxis += thisNucleus->semiaxes[j];
						}
						shortest_edist_tensor = userInputs.nucleation_parameters_list[userInputs.nucleation_parameters_list_index.at(thisNucleus->orderParameterIndex)].rotation_matrix * shortest_edist_tensor;
						for (unsigned int j=0; j<dim; j++){
							shortest_edist_tensor[j] /= userInputs.nucleation_parameters_list[userInputs.nucleation_parameters_list_index.at(thisNucleus->orderParameterIndex)].semiaxes[j];
						}
						r = shortest_edist_tensor.norm();
						avg_semiaxis /= dim;

						if (userInputs.nucleation_parameters_list_index.at(thisNucleus->orderParameterIndex) == 1){
							source_term1[i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/interface_coeff));
						}
						else {
							source_term2[i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/interface_coeff));
						}
					}
				}
			}
		}
	}
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
