
// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){
	// Variable 0
	// set_variable_name				(0,"f_tot");
	// set_variable_type				(0,SCALAR);
	//
	// set_need_value_residual_term	(0,true);
	// set_need_gradient_residual_term	(0,false);
	//
    // set_output_integral         	(0,true);

}

// =================================================================================
// Define the expressions for the post-processed fields
// =================================================================================

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//u
vectorgradType ux = variable_list.get_vector_gradient(0);

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
computeStress<dim>(CIJ, E, S);

scalarvalueType f_tot = constV(0.0);

for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  f_tot += constV(0.5) * S[i][j]*E[i][j];
  }
}


// Residuals for the equation to evolve the order parameter 
pp_variable_list.set_scalar_value_residual_term(0, f_tot);

}
