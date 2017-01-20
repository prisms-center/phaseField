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

  //Threads::ThreadMutex assembler_lock;
  //assembler_lock.acquire ();
  //energy+=cellEnergy;
  //assembler_lock.release ();
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



