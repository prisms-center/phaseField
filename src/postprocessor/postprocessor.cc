#include "../../include/postprocessor.h"
#include "../../include/list_of_CIJ.h"

template <int dim,int degree>
PostProcessor<dim,degree>::PostProcessor(userInputParameters<dim> _userInputs):userInputs(_userInputs){
}

template <int dim,int degree>
void PostProcessor<dim,degree>::computePostProcessedFields(MatrixFree<dim,double> matrixFreeObject, const std::vector<vectorType*> &solutionSet, std::vector<vectorType*> &postProcessedSet){


	// Zero out the postProcessedSet
	for(unsigned int fieldIndex=0; fieldIndex<userInputs.pp_number_of_variables; fieldIndex++){
		vectorType *U;
		U=new vectorType;
		postProcessedSet.push_back(U);
		matrixFreeObject.initialize_dof_vector(*U,  0); *U=0;
	}

	  //call to integrate and assemble
	  matrixFreeObject.cell_loop (&PostProcessor<dim,degree>::getPostProcessedFields, this, postProcessedSet, solutionSet);

}

template <int dim,int degree>
void PostProcessor<dim,degree>::getPostProcessedFields(const MatrixFree<dim,double> &data,
		std::vector<vectorType*> &dst,
		const std::vector<vectorType*> &src,
		const std::pair<unsigned int,unsigned int> &cell_range){

	//initialize FEEvaulation objects
	std::vector<dealii::FEEvaluation<dim,degree,degree+1,1,double> > scalar_vars;
	std::vector<dealii::FEEvaluation<dim,degree,degree+1,dim,double> > vector_vars;
	std::vector<dealii::FEEvaluation<dim,degree,degree+1,1,double> > pp_scalar_vars;
	std::vector<dealii::FEEvaluation<dim,degree,degree+1,dim,double> > pp_vector_vars;

	for (unsigned int i=0; i<userInputs.pp_number_of_variables; i++){
		if (userInputs.pp_varInfoList[i].is_scalar){
			dealii::FEEvaluation<dim,degree,degree+1,1,double> var(data, 0);
			pp_scalar_vars.push_back(var);
		}
		else {
			dealii::FEEvaluation<dim,degree,degree+1,dim,double> var(data, 0);
			pp_vector_vars.push_back(var);
		}
	}

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
	modelResidualsList.reserve(userInputs.pp_number_of_variables);

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

		// Initialize, read DOFs, and set evaulation flags for each variable
		for (unsigned int i=0; i<userInputs.pp_number_of_variables; i++){
			if (userInputs.pp_varInfoList[i].is_scalar) {
				pp_scalar_vars[userInputs.pp_varInfoList[i].scalar_or_vector_index].reinit(cell);
			}
			else {
				pp_vector_vars[userInputs.pp_varInfoList[i].scalar_or_vector_index].reinit(cell);
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
			postProcessedFields<dim>(modelVarList,modelResidualsList,q_point_loc,this->userInputs.material_moduli);

			// Submit values
			for (unsigned int i=0; i<userInputs.pp_number_of_variables; i++){
				if (userInputs.pp_varInfoList[i].is_scalar) {
					if (userInputs.pp_value_residual[i] == true){
						pp_scalar_vars[userInputs.pp_varInfoList[i].scalar_or_vector_index].submit_value(modelResidualsList[i].scalarValueResidual,q);
					}
					if (userInputs.pp_gradient_residual[i] == true){
						pp_scalar_vars[userInputs.pp_varInfoList[i].scalar_or_vector_index].submit_gradient(modelResidualsList[i].scalarGradResidual,q);
					}
				}
				else {
					if (userInputs.pp_value_residual[i] == true){
						pp_vector_vars[userInputs.pp_varInfoList[i].scalar_or_vector_index].submit_value(modelResidualsList[i].vectorValueResidual,q);
					}
					if (userInputs.pp_gradient_residual[i] == true){
						pp_vector_vars[userInputs.pp_varInfoList[i].scalar_or_vector_index].submit_gradient(modelResidualsList[i].vectorGradResidual,q);
					}
				}
			}

		}

		for (unsigned int i=0; i<userInputs.pp_number_of_variables; i++){
			if (userInputs.pp_varInfoList[i].is_scalar) {
				pp_scalar_vars[userInputs.pp_varInfoList[i].scalar_or_vector_index].integrate(userInputs.pp_value_residual[i], userInputs.pp_gradient_residual[i]);
				pp_scalar_vars[userInputs.pp_varInfoList[i].scalar_or_vector_index].distribute_local_to_global(*dst[userInputs.pp_varInfoList[i].global_var_index]);
			}
			else {
				pp_vector_vars[userInputs.pp_varInfoList[i].scalar_or_vector_index].integrate(userInputs.pp_value_residual[i], userInputs.pp_gradient_residual[i]);
				pp_vector_vars[userInputs.pp_varInfoList[i].scalar_or_vector_index].distribute_local_to_global(*dst[userInputs.pp_varInfoList[i].global_var_index]);
			}
		}
	}
}



template class PostProcessor<2,1>;
template class PostProcessor<3,1>;
template class PostProcessor<2,2>;
template class PostProcessor<3,2>;
template class PostProcessor<2,3>;
template class PostProcessor<3,3>;
