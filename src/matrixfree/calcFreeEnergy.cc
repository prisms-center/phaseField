//computeEnergy() method for MatrixFreePDE class

#ifndef COMPUTEENERGY_MATRIXFREE_H
#define COMPUTEENERGY_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//update RHS of each field
template <int dim>
void MatrixFreePDE<dim>::computeEnergy(){
  //log time
  computing_timer.enter_section("matrixFreePDE: computeEnergy");

  //call to integrate and assemble
  energy=0.0;
  energy_components.clear();
  energy_components.push_back(0.0);
  energy_components.push_back(0.0);
  energy_components.push_back(0.0);

  matrixFreeObject.cell_loop (&MatrixFreePDE<dim>::getEnergy, this, residualSet, solutionSet);

  //add across all processors
  energy=Utilities::MPI::sum(energy, MPI_COMM_WORLD);
  for (unsigned int i=0; i<3; i++){
	  energy_components[i]=Utilities::MPI::sum(energy_components[i], MPI_COMM_WORLD);
  }
  pcout << "Energy: " << energy << std::endl;
  pcout << "Energy Components: " << energy_components[0] << " " << energy_components[1] << " " << energy_components[2] << " " << std::endl;
  freeEnergyValues.push_back(energy);
  //end log
  computing_timer.exit_section("matrixFreePDE: computeEnergy");
}

template <int dim>
void  MatrixFreePDE<dim>::getEnergy(const MatrixFree<dim,double> &data,
				    std::vector<vectorType*> &dst,
				    const std::vector<vectorType*> &src,
				    const std::pair<unsigned int,unsigned int> &cell_range) {

	//initialize FEEvaulation objects
		  std::vector<typeScalar> scalar_vars;
		  std::vector<typeVector> vector_vars;

		  for (unsigned int i=0; i<num_var; i++){
			  if (varInfoListRHS[i].is_scalar){
				  typeScalar var(data, i);
				  scalar_vars.push_back(var);
			  }
			  else {
				  typeVector var(data, i);
				  vector_vars.push_back(var);
			  }
		  }

		  std::vector<modelVariable<dim> > modelVarList;
		  std::vector<modelResidual<dim> > modelResidualsList;
		  modelVarList.reserve(num_var);
		  modelResidualsList.reserve(num_var);

		  //loop over cells
		  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){

			  // Initialize, read DOFs, and set evaulation flags for each variable
			  for (unsigned int i=0; i<num_var; i++){
				  if (varInfoListRHS[i].is_scalar) {
					  scalar_vars[varInfoListRHS[i].scalar_or_vector_index].reinit(cell);
					  scalar_vars[varInfoListRHS[i].scalar_or_vector_index].read_dof_values_plain(*src[varInfoListRHS[i].global_var_index]);
					  scalar_vars[varInfoListRHS[i].scalar_or_vector_index].evaluate(need_value[i], need_gradient[i], need_hessian[i]);
				  }
				  else {
					  vector_vars[varInfoListRHS[i].scalar_or_vector_index].reinit(cell);
					  vector_vars[varInfoListRHS[i].scalar_or_vector_index].read_dof_values_plain(*src[varInfoListRHS[i].global_var_index]);
					  vector_vars[varInfoListRHS[i].scalar_or_vector_index].evaluate(need_value[i], need_gradient[i], need_hessian[i]);
				  }
			  }


			  unsigned int num_q_points;
			  if (scalar_vars.size() > 0){
				  num_q_points = scalar_vars[0].n_q_points;
			  }
			  else {
				  num_q_points = vector_vars[0].n_q_points;
			  }

			  dealii::AlignedVector<dealii::VectorizedArray<double> > JxW(num_q_points);

			  if (scalar_vars.size() > 0){
				  scalar_vars[0].fill_JxW_values(JxW);
			  }
			  else {
				  vector_vars[0].fill_JxW_values(JxW);
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

				  for (unsigned int i=0; i<num_var; i++){
					  if (varInfoListRHS[i].is_scalar) {
						  if (need_value[i]){
							  modelVarList[i].scalarValue = scalar_vars[varInfoListRHS[i].scalar_or_vector_index].get_value(q);
						  }
						  if (need_gradient[i]){
							  modelVarList[i].scalarGrad = scalar_vars[varInfoListRHS[i].scalar_or_vector_index].get_gradient(q);
						  }
						  if (need_hessian[i]){
							  modelVarList[i].scalarHess = scalar_vars[varInfoListRHS[i].scalar_or_vector_index].get_hessian(q);
						  }
					  }
					  else {
						  if (need_value[i]){
							  modelVarList[i].vectorValue = vector_vars[varInfoListRHS[i].scalar_or_vector_index].get_value(q);
						  }
						  if (need_gradient[i]){
							  modelVarList[i].vectorGrad = vector_vars[varInfoListRHS[i].scalar_or_vector_index].get_gradient(q);
						  }
						  if (need_hessian[i]){
							  modelVarList[i].vectorHess = vector_vars[varInfoListRHS[i].scalar_or_vector_index].get_hessian(q);
						  }
					  }
				  }

				  // Calculate the energy density
				  energyDensity(modelVarList,JxW[q],q_point_loc);
			  }
		  }

}

// output the integrated free energies into a text file
template <int dim>
void MatrixFreePDE<dim>::outputFreeEnergy(std::vector<double>& freeEnergyValues){

	  std::ofstream output_file("./freeEnergy.txt");
	  output_file.precision(10);
	  std::ostream_iterator<double> output_iterator(output_file, "\n");
	  std::copy(freeEnergyValues.begin(), freeEnergyValues.end(), output_iterator);
}




#endif



