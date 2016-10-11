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
					    *(ConstraintMatrix*) this->constraintsDirichletSet[currentFieldIndex]);
}

// Based on the contents of BC_list, mark faces on the triangulation as periodic
template <int dim>
void MatrixFreePDE<dim>::setPeriodicity(){
	// Default null implementation
}

// Set constraints to enforce periodic boundary conditions
template <int dim>
void MatrixFreePDE<dim>::setPeriodicityConstraints(ConstraintMatrix * constraints, DoFHandler<dim>* dof_handler){
	// Default null implementation
}

// Set constraints to pin the solution if there are no Dirichlet BCs for a component of a variable
template <int dim>
void MatrixFreePDE<dim>::setTranslationPreventionConstraints(ConstraintMatrix * constraints, DoFHandler<dim>* dof_handler){
	// Default null implementation
}

#endif
