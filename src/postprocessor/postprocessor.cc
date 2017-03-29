#include "../../include/postprocessor.h"

template <int dim,int degree>
PostProcessor<dim,degree>::PostProcessor(userInputParameters _userInputs):userInputs(_userInputs){
}

template <int dim,int degree>
void PostProcessor<dim,degree>::computePostProcessedFields(MatrixFree<dim,double> matrixFreeObject, const std::vector<vectorType*> &solutionSet, std::vector<vectorType*> &postProcessedSet){


	// Zero out the postProcessedSet
	for(unsigned int fieldIndex=0; fieldIndex<solutionSet.size(); fieldIndex++){
		vectorType *U;
		U=new vectorType;
		postProcessedSet.push_back(U);
		matrixFreeObject.initialize_dof_vector(*U,  fieldIndex); *U=0;


	}

	  //call to integrate and assemble
	  matrixFreeObject.cell_loop (&PostProcessor<dim,degree>::getPostProcessedFields, this, postProcessedSet, solutionSet);

	  std::cout << (*postProcessedSet[0])[10] << std::endl;

	  for(unsigned int fieldIndex=0; fieldIndex<solutionSet.size(); fieldIndex++){
		  postProcessedSet[fieldIndex]->update_ghost_values();
	  }

	  std::cout << (*postProcessedSet[0])[10] << std::endl;
}

template <int dim,int degree>
void PostProcessor<dim,degree>::getPostProcessedFields(const MatrixFree<dim,double> &data,
		std::vector<vectorType*> &dst,
		const std::vector<vectorType*> &src,
		const std::pair<unsigned int,unsigned int> &cell_range){

	//initialize FEEvaulation objects
	std::vector<dealii::FEEvaluation<dim,degree,degree+1,1,double>> scalar_vars;
	std::vector<dealii::FEEvaluation<dim,degree,degree+1,dim,double>> vector_vars;

	for (unsigned int i=0; i<userInputs.number_of_variables; i++){
		if (userInputs.varInfoListRHS[i].is_scalar){
			dealii::FEEvaluation<dim,degree,degree+1,1,double> var(data, i);
			scalar_vars.push_back(var);
		}
		else {
			dealii::FEEvaluation<dim,degree,degree+1,dim,double> var(data, i);
			vector_vars.push_back(var);
		}
	}

	std::vector<modelVariable<dim> > modelVarList;
	std::vector<modelResidual<dim> > modelResidualsList;
	modelVarList.reserve(userInputs.number_of_variables);
	modelResidualsList.reserve(userInputs.number_of_variables);

	//loop over cells
	for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){

		// Initialize, read DOFs, and set evaulation flags for each variable
		for (unsigned int i=0; i<userInputs.number_of_variables; i++){
			if (userInputs.varInfoListRHS[i].is_scalar) {
				scalar_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].reinit(cell);
				scalar_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].read_dof_values_plain(*src[userInputs.varInfoListRHS[i].global_var_index]);
				scalar_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].evaluate(userInputs.need_value[i], userInputs.need_gradient[i], userInputs.need_hessian[i]);
			}
			else {
				vector_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].reinit(cell);
				vector_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].read_dof_values_plain(*src[userInputs.varInfoListRHS[i].global_var_index]);
				vector_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].evaluate(userInputs.need_value[i], userInputs.need_gradient[i], userInputs.need_hessian[i]);
			}
		}

		unsigned int num_q_points;
		if (scalar_vars.size() > 0){
			num_q_points = scalar_vars[0].n_q_points;
		}
		else {
			num_q_points = vector_vars[0].n_q_points;
		}

		//loop over quadrature points
		for (unsigned int q=0; q<num_q_points; ++q){

			dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc;
			if (scalar_vars.size() > 0){
				q_point_loc = scalar_vars[0].quadrature_point(q);
			}
			else {
				q_point_loc = vector_vars[0].quadrature_point(q);
			}

			for (unsigned int i=0; i<userInputs.number_of_variables; i++){
				if (userInputs.varInfoListRHS[i].is_scalar) {
					if (userInputs.need_value[i]){
						modelVarList[i].scalarValue = scalar_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].get_value(q);
					}
					if (userInputs.need_gradient[i]){
						modelVarList[i].scalarGrad = scalar_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].get_gradient(q);
					}
					if (userInputs.need_hessian[i]){
						modelVarList[i].scalarHess = scalar_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].get_hessian(q);
					}
				}
				else {
					if (userInputs.need_value[i]){
						modelVarList[i].vectorValue = vector_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].get_value(q);
					}
					if (userInputs.need_gradient[i]){
						modelVarList[i].vectorGrad = vector_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].get_gradient(q);
					}
					if (userInputs.need_hessian[i]){
						modelVarList[i].vectorHess = vector_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].get_hessian(q);
					}
				}
			}

			// Calculate the residuals
			postProcessedFields(modelVarList,modelResidualsList,q_point_loc);

			// Submit values
			for (unsigned int i=0; i<userInputs.number_of_variables; i++){
				if (userInputs.varInfoListRHS[i].is_scalar) {
					if (userInputs.value_residual[i] == true){
						scalar_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].submit_value(modelResidualsList[i].scalarValueResidual,q);
					}
					if (userInputs.gradient_residual[i] == true){
						scalar_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].submit_gradient(modelResidualsList[i].scalarGradResidual,q);
					}
				}
				else {
					if (userInputs.value_residual[i] == true){
						vector_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].submit_value(modelResidualsList[i].vectorValueResidual,q);
					}
					if (userInputs.gradient_residual[i] == true){
						vector_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].submit_gradient(modelResidualsList[i].vectorGradResidual,q);
					}
				}
			}

		}

		for (unsigned int i=0; i<userInputs.number_of_variables; i++){
			if (userInputs.varInfoListRHS[i].is_scalar) {
				scalar_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].integrate(userInputs.value_residual[i], userInputs.gradient_residual[i]);
				scalar_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].distribute_local_to_global(*dst[userInputs.varInfoListRHS[i].global_var_index]);
			}
			else {
				vector_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].integrate(userInputs.value_residual[i], userInputs.gradient_residual[i]);
				vector_vars[userInputs.varInfoListRHS[i].scalar_or_vector_index].distribute_local_to_global(*dst[userInputs.varInfoListRHS[i].global_var_index]);
			}
		}
	}
}

template <int dim,int degree>
void PostProcessor<dim,degree>::postProcessedFields(const std::vector<modelVariable<dim> > & modelVariablesList,
												std::vector<modelResidual<dim> > & modelResidualsList,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// The order parameter and its derivatives (names here should match those in the macros above)
dealii::VectorizedArray<double> n = modelVariablesList[0].scalarValue;
dealii::Tensor<1, dim, dealii::VectorizedArray<double> > nx = modelVariablesList[0].scalarGrad;


// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
modelResidualsList[0].scalarValueResidual = 2.0*n; //constV(0.0);
modelResidualsList[0].scalarGradResidual = constV(0.0)*nx;

}

template class PostProcessor<2,1>;
template class PostProcessor<3,1>;
template class PostProcessor<2,2>;
template class PostProcessor<3,2>;
