//methods to apply boundary conditons 

#ifndef BOUNDARYCONDITIONS_MATRIXFREE_H
#define BOUNDARYCONDITIONS_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//methods to apply dirichlet BC's
template <int dim>
void MatrixFreePDE<dim>::applyDirichletBCs(){
  //default method to apply zero Dirichlet BC's on all components
  //of this field, given by currentFieldIndex, on all boundary faces
  VectorTools::interpolate_boundary_values (*this->dofHandlersSet[currentFieldIndex], \
					    0, ZeroFunction<dim>(fields[currentFieldIndex].numComponents), \
					    *(ConstraintMatrix*) this->constraintsSet[currentFieldIndex]);
}

#endif
