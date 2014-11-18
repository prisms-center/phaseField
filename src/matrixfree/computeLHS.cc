//computeLHS() method for MatrixFreePDE class temporarily treating

#ifndef COMPUTELHS_MATRIXFREE_H
#define COMPUTELHS_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized


//vmult operation for LHS
template <int dim>
void MatrixFreePDE<dim>::vmult (vectorType &dst, const vectorType &src) const{
  dst=0.0;
  matrixFreeObject.cell_loop (&MatrixFreePDE::computeLHS, this, dst, src);
  
  //Account for dirichlet BC's
  /*
  const std::vector<unsigned int>& constrained_dofs = data.get_constrained_dofs();
  for (unsigned int i=0; i<constrained_dofs.size(); ++i){
    unsigned int index=data.get_vector_partitioner()->local_to_global(constrained_dofs[i]);
    dst(index) += src(index);
  }
  */
}


//compute LHS
template <int dim>
void  MatrixFreePDE<dim>::computeLHS(const MatrixFree<dim,double> &data, 
				     vectorType &dst, 
				     const vectorType &src,
				     const std::pair<unsigned int,unsigned int> &cell_range) const{
  //initialize vals vectors
  unsigned int fieldIndex=0;
  if (fields[fieldIndex].type==SCALAR){
    typeScalar  vals(data,fieldIndex);
    
    //loop over all "cells"
    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
      //read values from corresponding solution vectors
      vals.reinit(cell);
      vals.read_dof_values_plain(src);
      vals.evaluate( getValue.find(fields[fieldIndex].name)->second, \
		     getGradient.find(fields[fieldIndex].name)->second, \
		     false);
    }
    //loop over quadrature points
    for (unsigned int q=0; q<vals.n_q_points; ++q){
      //getLHS(vals, q);
    }
    vals.integrate( setValue.find(fields[fieldIndex].name)->second,	\
		    setGradient.find(fields[fieldIndex].name)->second);
    vals.distribute_local_to_global(dst); 
  }
  else{
    typeVector  vals(data,fieldIndex);
    
    //loop over all "cells"
    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
      //read values from corresponding solution vectors
      vals.reinit(cell);
      vals.read_dof_values_plain(src);
      vals.evaluate( getValue.find(fields[fieldIndex].name)->second,	\
		     getGradient.find(fields[fieldIndex].name)->second, \
		     false);
    }
    //loop over quadrature points
    for (unsigned int q=0; q<vals.n_q_points; ++q){
      //getLHS(vals, q);
    }
    vals.integrate( setValue.find(fields[fieldIndex].name)->second,	\
		    setGradient.find(fields[fieldIndex].name)->second);
      vals.distribute_local_to_global(dst); 
  }
}

#endif

