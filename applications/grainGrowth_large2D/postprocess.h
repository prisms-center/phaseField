// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadPostProcessorVariableAttributes(){

	// Variable 0
	set_variable_name				(0,"feature_ids");
	set_variable_type				(0,SCALAR);

    set_dependencies_value_residual_term_RHS(0, "n0, n1, n2, n3, n4, n5, n6, n7");
    set_dependencies_gradient_residual_term_RHS(0, "");

    set_output_integral         	(0,true);

}

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

scalarvalueType ni;

scalarvalueType max_val = constV(-1.0);
dealii::VectorizedArray<int> max_op = constV(0);
for (unsigned int i=0; i<userInputs.number_of_variables; i++){
    ni = variable_list.get_scalar_value(i) ;

    for (unsigned int v=0; v<ni.n_array_elements;v++){
        if (ni[v] > max_val[v]){
            max_val[v] = ni[v];
            max_op[v] = i;
        }
    }
}

scalarvalueType feature_ids = constV(-1.0);
for (unsigned int v=0; v<ni.n_array_elements;v++){
    for (unsigned int g=0; g<this->simplified_grain_representations.size(); g++){
        if (this->simplified_grain_representations[g].getOldOrderParameterId() == max_op[v]){
            dealii::Point<dim> q_point_loc_nonvec;
            for (unsigned int d=0;d<dim;d++){
                q_point_loc_nonvec(d) = q_point_loc(d)[v];
            }

            if ( q_point_loc_nonvec.distance(this->simplified_grain_representations[g].getCenter()) < this->simplified_grain_representations[g].getRadius() ){
                feature_ids[v] = (double)(this->simplified_grain_representations[g].getGrainId());
            }
        }
    }
}

pp_variable_list.set_scalar_value_residual_term(0, feature_ids);

}
