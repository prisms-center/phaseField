#ifndef MODIFYFIELDS_MATRIXFREE_H
#define MODIFYFIELDS_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
void  MatrixFreePDE<dim,degree>::modifySolutionFields(){
  //default trivial implementation.
}

#ifndef MATRIXFREEPDE_TEMPLATE_INSTANTIATION
#define MATRIXFREEPDE_TEMPLATE_INSTANTIATION
template class MatrixFreePDE<2,1>;
template class MatrixFreePDE<3,1>;
#endif


#endif

