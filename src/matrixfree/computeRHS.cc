//computeRHS() method for MatrixFreePDE class

#ifndef COMPUTERHS_MATRIXFREE_H
#define COMPUTERHS_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//update RHS of each field
template <int dim>
void MatrixFreePDE<dim>::computeRHS(){
  //log time
  computing_timer.enter_section("matrixFreePDE: computeRHS");

  //clear residual vectors before update
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    (*residualSet[fieldIndex])=0.0;
  }

  //call to integrate and assemble 
  matrixFreeObject.cell_loop (&MatrixFreePDE<dim>::getRHS, this, residualSet, solutionSet);

  //end log
  computing_timer.exit_section("matrixFreePDE: computeRHS");
}

#endif

