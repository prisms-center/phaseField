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

  //create temporary copy of src vector as src2, as vector src is marked const and cannot be changed
  vectorType src2;
  matrixFreeObject.initialize_dof_vector(src2,  ellipticFieldIndex);
  src2=src;
  
  //set Dirichlet nodes force to zero in the src
  for (std::map<types::global_dof_index, double>::const_iterator it=valuesDirichletSet[currentFieldIndex]->begin(); it!=valuesDirichletSet[currentFieldIndex]->end(); ++it){
    if (src2.in_local_range(it->first)){
      src2(it->first) = 0.0; //*jacobianDiagonal(it->first);
    }
  }
  constraintsOtherSet[currentFieldIndex]->distribute(src2);

  //call cell_loop 
  dst=0.0;
  matrixFreeObject.cell_loop (&MatrixFreePDE<dim>::getLHS, this, dst, src2);
  dst.compress(VectorOperation::add);
  
  //Account for Dirichlet BC's (essentially copy dirichlet DOF values present in src to dst)
  for (std::map<types::global_dof_index, double>::const_iterator it=valuesDirichletSet[currentFieldIndex]->begin(); it!=valuesDirichletSet[currentFieldIndex]->end(); ++it){
    if (dst.in_local_range(it->first)){
      dst(it->first) = src(it->first); //*jacobianDiagonal(it->first);
    }
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


