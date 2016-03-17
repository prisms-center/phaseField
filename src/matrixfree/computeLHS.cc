//vmult() and getLHS() method for MatrixFreePDE class 
#ifndef COMPUTELHS_MATRIXFREE_H
#define COMPUTELHS_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//vmult operation for LHS
template <int dim>
void MatrixFreePDE<dim>::vmult (vectorType &dst, const vectorType &src) const{
  //log time
  computing_timer.enter_section("matrixFreePDE: computeLHS");
  
  dst=0.0;  
  matrixFreeObject.cell_loop (&MatrixFreePDE<dim>::getLHS, this, dst, src);
  
  //Account for dirichlet BC's (essentially copy dirichlet DOF values present in src to dst)
  const std::vector<unsigned int>& constrained_dofs = matrixFreeObject.get_constrained_dofs(currentFieldIndex);
  for (unsigned int i=0; i<constrained_dofs.size(); ++i){
    unsigned int globalIndex = matrixFreeObject.get_vector_partitioner(currentFieldIndex)->local_to_global(constrained_dofs[i]);
    dst(globalIndex) += src(globalIndex); //note: check if "+" required
  } 
  
  //end log
  computing_timer.exit_section("matrixFreePDE: computeLHS");
}
  
template <int dim>
void  MatrixFreePDE<dim>::getLHS(const MatrixFree<dim,double> &data, 
				 vectorType &dst, 
				 const vectorType &src,
				 const std::pair<unsigned int,unsigned int> &cell_range) const{
  pcout << "\n\nError: computeLHS.cc: getLHS(typeScalar& , unsigned int ) not implemented in the derived class, but is called\n";
  exit(-1);
}

#endif


