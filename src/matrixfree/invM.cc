//computeInvM() method for MatrixFreePDE class

#ifndef INVM_MATRIXFREE_H
#define INVM_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//compute inverse of the diagonal mass matrix and store in vector invM
template <int dim>
void MatrixFreePDE<dim>::computeInvM(){
  //initialize  invM
  bool invMInitialized=false;
  unsigned int parabolicFieldIndex=0;
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    if (fields[fieldIndex].pdetype==PARABOLIC){
      matrixFreeObject.initialize_dof_vector (invM, fieldIndex);
      parabolicFieldIndex=fieldIndex;
      invMInitialized=true;
      break;
    }
  }
  //check if invM initialized
  if (!invMInitialized){
    pcout << "matrixFreePDE.h: no PARABOLIC field... hence setting parabolicFieldIndex to 0 and marching ahead withn invM computation\n";
    //exit(-1);
  }
  
  //compute invM
  matrixFreeObject.initialize_dof_vector (invM, parabolicFieldIndex);
  invM=0.0;
  VectorizedArray<double> one = make_vectorized_array (1.0);
  
  //select gauss lobatto quad points which are suboptimal but give diogonal M 
  FEEvaluation<dim,finiteElementDegree> fe_eval(matrixFreeObject, parabolicFieldIndex);
  const unsigned int n_q_points = fe_eval.n_q_points;
  for (unsigned int cell=0; cell<matrixFreeObject.n_macro_cells(); ++cell){
    fe_eval.reinit(cell);
    for (unsigned int q=0; q<n_q_points; ++q){
      fe_eval.submit_value(one,q);
    }
    fe_eval.integrate (true,false);
    fe_eval.distribute_local_to_global (invM);
  }
  invM.compress(VectorOperation::add);
  
  //invert mass matrix diagonal elements
  for (unsigned int k=0; k<invM.local_size(); ++k){
    if (std::abs(invM.local_element(k))>1.0e-15){
      invM.local_element(k) = 1./invM.local_element(k);
    }
    else{
      invM.local_element(k) = 0;
    }
  } 
  pcout << "computed mass matrix (using FE space for field: " << parabolicFieldIndex << ")\n";
}

#endif
