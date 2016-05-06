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
  matrixFreeObject.cell_loop (&MatrixFreePDE<dim>::getEnergy, this, residualSet, solutionSet);
  //add across all processors
  energy=Utilities::MPI::sum(energy, MPI_COMM_WORLD);
  pcout << "Energy: " << energy << std::endl;
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


#endif



