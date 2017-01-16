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

// Determine which (if any) components of the current field have rigid body modes (i.e no Dirichlet BCs) if the
// equation is elliptic
template <int dim>
void MatrixFreePDE<dim>::getComponentsWithRigidBodyModes( std::vector<int> & rigidBodyModeComponents){

}

//// Set constraints to pin the solution if there are no Dirichlet BCs for a component of a variable
//template <int dim>
//void MatrixFreePDE<dim>::setRigidBodyModeConstraints(std::vector<int> rigidBodyModeComponents, ConstraintMatrix * constraints, DoFHandler<dim>* dof_handler){
//	// Default null implementation
//}

// Set constraints to pin the solution if there are no Dirichlet BCs for a component of a variable in an elliptic equation
template <int dim>
void MatrixFreePDE<dim>::setRigidBodyModeConstraints( std::vector<int> rigidBodyModeComponents, ConstraintMatrix * constraints, DoFHandler<dim>* dof_handler){

	if ( rigidBodyModeComponents.size() > 0 ){

		// Choose the point where the constraint will be placed. Must be the coordinates of a vertex.
		dealii::Point<dim> target_point(0,0);

		unsigned int vertices_per_cell=GeometryInfo<dim>::vertices_per_cell;

		// Loop over each locally owned cell
		typename DoFHandler<dim>::active_cell_iterator cell= dof_handler->begin_active(), endc = dof_handler->end();

		for (; cell!=endc; ++cell){
			if (cell->is_locally_owned()){
				for (unsigned int i=0; i<vertices_per_cell; ++i){

					// Check if the vertex is the target vertex
					if (target_point.distance (cell->vertex(i)) < 1e-2 * cell->diameter()){

						// Loop through the list of components with rigid body modes and add an inhomogeneous constraint for each
						for (unsigned int component_num = 0; component_num < rigidBodyModeComponents.size(); component_num++){
							unsigned int nodeID=cell->vertex_dof_index(i,component_num);
							constraints->add_line(nodeID);
							constraints->set_inhomogeneity(nodeID,0.0);
						}
				   }
			   }
		   }
	   }
   }
}

#endif
