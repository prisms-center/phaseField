//methods to mark boundaries 

#ifndef MARKBOUNDARIES_MATRIXFREE_H
#define MARKBOUNDARIES_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//methods to mark boundaries
template <int dim>
void MatrixFreePDE<dim>::markBoundaries(){
  //default method does nothing, as by default all boundary faces are
  //marked with flag 0
}

#endif
