//computeLHS() method for MatrixFreePDE class temporarily treating

#ifndef COMPUTELHS_MATRIXFREE_H
#define COMPUTELHS_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//vmult operation for LHS
template <int dim>
void MatrixFreePDE<dim>::vmult (vectorType &dst, const vectorType &src) const{
  //log time
  computing_timer.enter_section("matrixFreePDE: updateLHS (vmult)");
  
  dst=0.0;  
  //scalar field
  if (fields[currentFieldIndex].type==SCALAR){
    matrixFreeObject.cell_loop (&MatrixFreePDE<dim>::computeLHS<typeScalar>, this, dst, src);
  }
  //vector field
  else if (fields[currentFieldIndex].type==VECTOR){
    matrixFreeObject.cell_loop (&MatrixFreePDE<dim>::computeLHS<typeVector>, this, dst, src);
  }
  
  //Account for dirichlet BC's (essentially copy dirichlet DOF values present in src to dst)
  const std::vector<unsigned int>& constrained_dofs = matrixFreeObject.get_constrained_dofs(currentFieldIndex);
  for (unsigned int i=0; i<constrained_dofs.size(); ++i){
    unsigned int globalIndex = matrixFreeObject.get_vector_partitioner(currentFieldIndex)->local_to_global(constrained_dofs[i]);
    dst(globalIndex) += src(globalIndex); //note: check if "+" required
  } 
  
  //end log
  computing_timer.exit_section("matrixFreePDE: updateLHS (vmult)");
}

//compute LHS.  addtional template typename "T" is to identify whether
//its being called for a scalar or vector field, i.e, to choose
//typeScalar or typeVector for vals
template <int dim>
template <typename T>
void  MatrixFreePDE<dim>::computeLHS(const MatrixFree<dim,double> &data, 
				     vectorType &dst, 
				     const vectorType &src,
				     const std::pair<unsigned int,unsigned int> &cell_range) const{
  //initialize vals vectors
  T  vals(data,currentFieldIndex);
  Assert((vals.n_q_points) == (num_quadrature_points), ::ExcDimensionMismatch((vals.n_q_points),(num_quadrature_points))); //replace with a better Exception
  
  //loop over all "cells" (all super-cells) (note: cell_range.second is the last super-cell allocated to this thread)
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //read values from corresponding solution vectors
    vals.reinit(cell);
    vals.read_dof_values(src);
    vals.evaluate( getValue.find(fields[currentFieldIndex].name)->second, \
		   getGradient.find(fields[currentFieldIndex].name)->second, \
		   false);
  
    //loop over quadrature points
    for (unsigned int q=0; q<num_quadrature_points; ++q){
      getLHS(vals, q);
    }
    
    //integrate
    vals.integrate( setValue.find(fields[currentFieldIndex].name)->second, \
		    setGradient.find(fields[currentFieldIndex].name)->second);
    
    //assemble
    vals.distribute_local_to_global(dst); 
  }
}

template <int dim>
void  MatrixFreePDE<dim>::getLHS(typeScalar& vals, unsigned int q) const{
  pcout << "\ncomputeLHS.cc: getLHS(typeScalar& , unsigned int ) not implemented in the derived class, but is called\n";
  exit(-1);
}

template <int dim>
void  MatrixFreePDE<dim>::getLHS(typeVector& vals, unsigned int q) const{
  pcout << "\ncomputeLHS.cc: getLHS(typeVector& , unsigned int ) not implemented in the derived class, but is called\n";
  exit(-1);
}


#endif

