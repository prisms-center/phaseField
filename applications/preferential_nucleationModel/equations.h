
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
	set_variable_name				(1,"n");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,PARABOLIC);

	set_need_value					(1,true);
	set_need_gradient				(1,true);
	set_need_hessian				(1,false);

	set_need_value_residual_term	(1,true);
	set_need_gradient_residual_term	(1,true);

	set_allowed_to_nucleate			(1, true);
	set_need_value_nucleation		(1, true);
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
#define hV (3.0*n*n - 2.0*n*n*n)
#define hnV (6.0*n - 6.0*n*n)

// KKS model c_alpha and c_beta as a function of c and h
#define c_alpha ((B2*(c-cbtmin*hV) + A2*calmin*hV)/(A2*hV+B2*(1.0-hV)))
#define c_beta ((A2*(c-calmin*(1.0-hV))+B2*cbtmin*(1.0-hV))/(A2*hV+B2*(1.0-hV)))

// Double-Well function (can be used to tune the interfacial energy)
#define fbarrierV (n*n - 2.0*n*n*n + n*n*n*n)
#define fbarriernV (2.0*n - 6.0*n*n + 4.0*n*n*n)

// Residual equations
// For concentration
#define term_muxV (cx + (c_alpha - c_beta)*hnV*nx)
#define rcV   (c)
#define rcxV  (constV(-McV*userInputs.dtValue)*term_muxV)
//For order parameter (gamma is a variable order parameter mobility factor)
#define rnV   (n-constV(userInputs.dtValue*MnV)*gamma*((fbV-faV)*hnV - (c_beta-c_alpha)*fbcV*hnV + W_barrier*fbarriernV))
#define rnxV  (constV(-userInputs.dtValue*KnV*MnV)*gamma*nx)

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
scalarvalueType n = variable_list.get_scalar_value(1);
scalargradType nx = variable_list.get_scalar_gradient(1);

// -------------------------------------------------
// Nucleation expressions
// -------------------------------------------------
dealii::VectorizedArray<double> source_term = constV(0.0);
dealii::VectorizedArray<double> gamma = constV(1.0);
seedNucleus(q_point_loc,source_term,gamma);
// -------------------------------------------------

// Residuals for the equation to evolve the concentration
variable_list.set_scalar_value_residual_term(0,rcV);
variable_list.set_scalar_gradient_residual_term(0,rcxV);

// Residuals for the equation to evolve the order parameter
variable_list.set_scalar_value_residual_term(1,rnV + source_term);
variable_list.set_scalar_gradient_residual_term(1,rnxV);

}

// =================================================================================
// seedNucleus: a function particular to this app
// =================================================================================
template <int dim,int degree>
void customPDE<dim,degree>::seedNucleus(const dealii::Point<dim, dealii::VectorizedArray<double> > & q_point_loc,
	dealii::VectorizedArray<double> & source_term,
	dealii::VectorizedArray<double> & gamma) const {

		// Loop through all of the seeded nuclei
		for (typename std::vector<nucleus<dim> >::const_iterator thisNucleus=this->nuclei.begin(); thisNucleus!=this->nuclei.end(); ++thisNucleus){

			if (thisNucleus->seededTime + thisNucleus->seedingTime > this->currentTime){

				// Calculate the weighted distance function to the order parameter freeze boundary (weighted_dist = 1.0 on that boundary)
				dealii::VectorizedArray<double> weighted_dist = this->weightedDistanceFromNucleusCenter(thisNucleus->center, userInputs.get_nucleus_freeze_semiaxes(thisNucleus->orderParameterIndex), q_point_loc, thisNucleus->orderParameterIndex);

				for (unsigned i=0; i<gamma.n_array_elements;i++){
					if (weighted_dist[i] <= 1.0){
						gamma[i] = 0.0;

						// Seed a nucleus if it was added to the list of nuclei this time step
						if (thisNucleus->seedingTimestep == this->currentIncrement){
							// Find the weighted distance to the outer edge of the nucleus and use it to calculate the order parameter source term (r = 1.0 on that boundary)
							dealii::Point<dim,double> q_point_loc_element;
							for (unsigned int j=0; j<dim; j++){
								q_point_loc_element(j) = q_point_loc(j)[i];
							}
							double r = this->weightedDistanceFromNucleusCenter(thisNucleus->center, userInputs.get_nucleus_semiaxes(thisNucleus->orderParameterIndex), q_point_loc_element, thisNucleus->orderParameterIndex);

							double avg_semiaxis = 0.0;
							for (unsigned int j=0; j<dim; j++){
								avg_semiaxis += thisNucleus->semiaxes[j];
							}
							avg_semiaxis /= dim;

							source_term[i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/interface_coeff));
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
