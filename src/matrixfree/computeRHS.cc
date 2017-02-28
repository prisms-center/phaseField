//computeRHS() method for MatrixFreePDE class

#ifndef COMPUTERHS_MATRIXFREE_H
#define COMPUTERHS_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../include/matrixFreePDE.h"

//update RHS of each field
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::computeRHS(){
  //log time
  computing_timer.enter_section("matrixFreePDE: computeRHS");

  //clear residual vectors before update
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    (*residualSet[fieldIndex])=0.0;
  }

  //call to integrate and assemble 
  matrixFreeObject.cell_loop (&MatrixFreePDE<dim,degree>::getRHS, this, residualSet, solutionSet);

  //end log
  computing_timer.exit_section("matrixFreePDE: computeRHS");
}

template <int dim, int degree>
void MatrixFreePDE<dim,degree>::getRHS(const MatrixFree<dim,double> &data,
					       std::vector<vectorType*> &dst,
					       const std::vector<vectorType*> &src,
					       const std::pair<unsigned int,unsigned int> &cell_range) const{


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
		  residualRHS(modelVarList,modelResidualsList,q_point_loc);

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

#ifndef MATRIXFREEPDE_TEMPLATE_INSTANTIATION
#define MATRIXFREEPDE_TEMPLATE_INSTANTIATION
template class MatrixFreePDE<2,1>;
template class MatrixFreePDE<3,1>;
#endif



#endif

