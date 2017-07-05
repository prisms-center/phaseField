//vmult() and getLHS() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

//vmult operation for LHS
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::vmult (vectorType &dst, const vectorType &src) const{
  //log time
  computing_timer.enter_section("matrixFreePDE: computeLHS");

  //create temporary copy of src vector as src2, as vector src is marked const and cannot be changed
  dealii::parallel::distributed::Vector<double> src2;
  matrixFreeObject.initialize_dof_vector(src2,  currentFieldIndex);
  src2=src;

  //call cell_loop
  dst=0.0;
  matrixFreeObject.cell_loop (&MatrixFreePDE<dim,degree>::getLHS, this, dst, src2);

  //Account for Dirichlet BC's (essentially copy dirichlet DOF values present in src to dst, although it is unclear why the constraints can't just be distributed here)
  for (std::map<types::global_dof_index, double>::const_iterator it=valuesDirichletSet[currentFieldIndex]->begin(); it!=valuesDirichletSet[currentFieldIndex]->end(); ++it){
    if (dst.in_local_range(it->first)){
      dst(it->first) = src(it->first); //*jacobianDiagonal(it->first);
    }
  }

  //end log
  computing_timer.exit_section("matrixFreePDE: computeLHS");
}

template <int dim, int degree>
void  MatrixFreePDE<dim,degree>::getLHS(const MatrixFree<dim,double> &data,
				 vectorType &dst,
				 const vectorType &src,
				 const std::pair<unsigned int,unsigned int> &cell_range) const{


	variable_info resInfoLHS;
	for (unsigned int i=0; i<userInputs.num_var_LHS; i++){
		if (currentFieldIndex == userInputs.varInfoListLHS[i].global_var_index){
			resInfoLHS = userInputs.varInfoListLHS[i];
			break;
		}
	}

	//initialize FEEvaulation objects
	std::vector<dealii::FEEvaluation<dim,degree,degree+1,1,double> > scalar_vars;
	std::vector<dealii::FEEvaluation<dim,degree,degree+1,dim,double> > vector_vars;

	for (unsigned int i=0; i<userInputs.num_var_LHS; i++){
		if (userInputs.varInfoListLHS[i].is_scalar){
			dealii::FEEvaluation<dim,degree,degree+1,1,double> var(data, userInputs.varInfoListLHS[i].global_var_index);
			scalar_vars.push_back(var);
		}
		else {
			dealii::FEEvaluation<dim,degree,degree+1,dim,double> var(data, userInputs.varInfoListLHS[i].global_var_index);
			vector_vars.push_back(var);
		}
	}

	std::vector<modelVariable<dim> > modelVarList;
	modelVarList.reserve(userInputs.num_var_LHS);
	modelResidual<dim> modelRes;

	//loop over cells
	for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){

		// Initialize, read DOFs, and set evaulation flags for each variable
		for (unsigned int i=0; i<userInputs.num_var_LHS; i++){
			if (userInputs.varInfoListLHS[i].is_scalar) {
				scalar_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].reinit(cell);
				if ( userInputs.varInfoListLHS[i].global_var_index == resInfoLHS.global_var_index ){
					scalar_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].read_dof_values(src);
				}
				else{
					scalar_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].read_dof_values(*solutionSet[userInputs.varInfoListLHS[i].global_var_index]);
				}
				scalar_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].evaluate(userInputs.need_value_LHS[userInputs.varInfoListLHS[i].global_var_index], userInputs.need_gradient_LHS[userInputs.varInfoListLHS[i].global_var_index], userInputs.need_hessian_LHS[userInputs.varInfoListLHS[i].global_var_index]);
			}
			else {
				vector_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].reinit(cell);
				if ( userInputs.varInfoListLHS[i].global_var_index == resInfoLHS.global_var_index ){
					vector_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].read_dof_values(src);
				}
				else {
					vector_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].read_dof_values(*solutionSet[userInputs.varInfoListLHS[i].global_var_index]);
				}
				vector_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].evaluate(userInputs.need_value_LHS[userInputs.varInfoListLHS[i].global_var_index], userInputs.need_gradient_LHS[userInputs.varInfoListLHS[i].global_var_index], userInputs.need_hessian_LHS[userInputs.varInfoListLHS[i].global_var_index]);
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

			for (unsigned int i=0; i<userInputs.num_var_LHS; i++){
				if (userInputs.varInfoListLHS[i].is_scalar) {
					if (userInputs.need_value_LHS[userInputs.varInfoListLHS[i].global_var_index]){
						modelVarList[i].scalarValue = scalar_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].get_value(q);
					}
					if (userInputs.need_gradient_LHS[userInputs.varInfoListLHS[i].global_var_index]){
						modelVarList[i].scalarGrad = scalar_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].get_gradient(q);
					}
					if (userInputs.need_hessian_LHS[userInputs.varInfoListLHS[i].global_var_index]){
						modelVarList[i].scalarHess = scalar_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].get_hessian(q);
					}
				}
				else {
					if (userInputs.need_value_LHS[userInputs.varInfoListLHS[i].global_var_index]){
						modelVarList[i].vectorValue = vector_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].get_value(q);
					}
					if (userInputs.need_gradient_LHS[userInputs.varInfoListLHS[i].global_var_index]){
						modelVarList[i].vectorGrad = vector_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].get_gradient(q);
					}
					if (userInputs.need_hessian_LHS[userInputs.varInfoListLHS[i].global_var_index]){
						modelVarList[i].vectorHess = vector_vars[userInputs.varInfoListLHS[i].scalar_or_vector_index].get_hessian(q);
					}
				}
			}

			// Calculate the residuals
			residualLHS(modelVarList,modelRes,q_point_loc);

			// Submit values
			if (resInfoLHS.is_scalar){
				if (userInputs.value_residual[resInfoLHS.global_var_index]){
					scalar_vars[resInfoLHS.scalar_or_vector_index].submit_value(modelRes.scalarValueResidual,q);
				}
				if (userInputs.gradient_residual[resInfoLHS.global_var_index]){
					scalar_vars[resInfoLHS.scalar_or_vector_index].submit_gradient(modelRes.scalarGradResidual,q);
				}
			}
			else {
				if (userInputs.value_residual[resInfoLHS.global_var_index]){
					vector_vars[resInfoLHS.scalar_or_vector_index].submit_value(modelRes.vectorValueResidual,q);
				}
				if (userInputs.gradient_residual[resInfoLHS.global_var_index]){
					vector_vars[resInfoLHS.scalar_or_vector_index].submit_gradient(modelRes.vectorGradResidual,q);
				}
			}

		}

		//integrate
		if (resInfoLHS.is_scalar) {
			scalar_vars[resInfoLHS.scalar_or_vector_index].integrate(userInputs.value_residual[resInfoLHS.global_var_index], userInputs.gradient_residual[resInfoLHS.global_var_index]);
			scalar_vars[resInfoLHS.scalar_or_vector_index].distribute_local_to_global(dst);
		}
		else {
			vector_vars[resInfoLHS.scalar_or_vector_index].integrate(userInputs.value_residual[resInfoLHS.global_var_index], userInputs.gradient_residual[resInfoLHS.global_var_index]);
			vector_vars[resInfoLHS.scalar_or_vector_index].distribute_local_to_global(dst);
		}
	}
}

#include "../../include/matrixFreePDE_template_instantiations.h"
