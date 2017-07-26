//methods to apply boundary conditons

//#ifndef's) till library packaging scheme is finalized

#include "../../include/matrixFreePDE.h"
#include "../../include/vectorBCFunction.h"

#include "../../include/nonUniformDirichletBC.h"


//methods to apply dirichlet BC's
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::applyDirichletBCs(){
	// First, get the variable index of the current field
	  unsigned int starting_BC_list_index = 0;

	  for (unsigned int i=0; i<currentFieldIndex; i++){

		  if (userInputs.var_type[i] == SCALAR){
			  starting_BC_list_index++;
		  }
		  else {
			  starting_BC_list_index+=dim;
		  }
	  }

	  if (userInputs.var_type[currentFieldIndex] == SCALAR){
		  for (unsigned int direction = 0; direction < 2*dim; direction++){
			  if (userInputs.BC_list[starting_BC_list_index].var_BC_type[direction] == DIRICHLET){
				  VectorTools::interpolate_boundary_values (*dofHandlersSet[currentFieldIndex],\
						  direction, ConstantFunction<dim>(userInputs.BC_list[starting_BC_list_index].var_BC_val[direction],1), *(ConstraintMatrix*) \
						  constraintsDirichletSet[currentFieldIndex]);

			  }
			  else if (userInputs.BC_list[starting_BC_list_index].var_BC_type[direction] == NON_UNIFORM_DIRICHLET){
				  VectorTools::interpolate_boundary_values (*dofHandlersSet[currentFieldIndex],\
  						direction, NonUniformDirichletBC<dim>(currentFieldIndex,direction,userInputs), *(ConstraintMatrix*) \
  						constraintsDirichletSet[currentFieldIndex]);
			  }
		  }
	  }
	  else {
		  for (unsigned int direction = 0; direction < 2*dim; direction++){

			  std::vector<double> BC_values;
			  for (unsigned int component=0; component < dim; component++){
				  BC_values.push_back(userInputs.BC_list[starting_BC_list_index+component].var_BC_val[direction]);
			  }

			  std::vector<bool> mask;
			  for (unsigned int component=0; component < dim; component++){
				  if (userInputs.BC_list[starting_BC_list_index+component].var_BC_type[direction] == DIRICHLET){
					  mask.push_back(true);
				  }
				  else {
					  mask.push_back(false);
				  }
			  }

			  VectorTools::interpolate_boundary_values (*dofHandlersSet[currentFieldIndex],\
					  direction, vectorBCFunction<dim>(BC_values), *(ConstraintMatrix*) \
					  constraintsDirichletSet[currentFieldIndex],mask);

				// Mask again, this time for non-uniform Dirichlet BCs
			  mask.clear();
			  for (unsigned int component=0; component < dim; component++){
				  if (userInputs.BC_list[starting_BC_list_index+component].var_BC_type[direction] == NON_UNIFORM_DIRICHLET){
					  mask.push_back(true);
				  }
				  else {
					  mask.push_back(false);
				  }
			  }

			  VectorTools::interpolate_boundary_values (*dofHandlersSet[currentFieldIndex],\
				  direction, NonUniformDirichletBCVec<dim>(currentFieldIndex,direction,userInputs), *(ConstraintMatrix*) \
				  constraintsDirichletSet[currentFieldIndex],mask);


		  }
	  }
}

// Based on the contents of BC_list, mark faces on the triangulation as periodic
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::setPeriodicity(){
	std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> > periodicity_vector;
		for (int i=0; i<dim; ++i){
			bool periodic_pair = false;
			for (unsigned int field_num=0; field_num < userInputs.BC_list.size(); field_num++){
				if (userInputs.BC_list[field_num].var_BC_type[2*i] == PERIODIC){
					periodic_pair = true;
				}
			}
			if (periodic_pair == true){
				GridTools::collect_periodic_faces(triangulation, /*b_id1*/ 2*i, /*b_id2*/ 2*i+1,
								/*direction*/ i, periodicity_vector);
			}
		}

		triangulation.add_periodicity(periodicity_vector);
		pcout << "periodic facepairs: " << periodicity_vector.size() << std::endl;
}

// Set constraints to enforce periodic boundary conditions
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::setPeriodicityConstraints(ConstraintMatrix * constraints, const DoFHandler<dim>* dof_handler) const {
	// First, get the variable index of the current field
		unsigned int starting_BC_list_index = 0;
		for (unsigned int i=0; i<currentFieldIndex; i++){
			if (userInputs.var_type[i] == SCALAR){
				starting_BC_list_index++;
			}
			else {
				starting_BC_list_index+=dim;
			}
		}

		std::vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator> > periodicity_vector;
	    for (int i=0; i<dim; ++i){
	    	if (userInputs.BC_list[starting_BC_list_index].var_BC_type[2*i] == PERIODIC){
	    		GridTools::collect_periodic_faces(*dof_handler, /*b_id1*/ 2*i, /*b_id2*/ 2*i+1,
	    				/*direction*/ i, periodicity_vector);
	    	}
	    }
	    DoFTools::make_periodicity_constraints<DoFHandler<dim> >(periodicity_vector, *constraints);
}

// Determine which (if any) components of the current field have rigid body modes (i.e no Dirichlet BCs) if the
// equation is elliptic
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::getComponentsWithRigidBodyModes( std::vector<int> & rigidBodyModeComponents) const {
	// Rigid body modes only matter for elliptic equations
		if (userInputs.var_eq_type[currentFieldIndex] == ELLIPTIC){

			// First, get the variable index of the current field
			unsigned int starting_BC_list_index = 0;
			for (unsigned int i=0; i<currentFieldIndex; i++){
				if (userInputs.var_type[i] == SCALAR){
					starting_BC_list_index++;
				}
				else {
					starting_BC_list_index+=dim;
				}
			}

			// Get number of components of the field
			unsigned int num_components = 1;
			if (userInputs.var_type[currentFieldIndex] == VECTOR){
				num_components = dim;
			}

			// Loop over each component and determine if it has a rigid body mode (i.e. no Dirichlet BCs)
			for (unsigned int component=0; component < num_components; component++){
				bool rigidBodyMode = true;
				for (unsigned int direction = 0; direction < 2*dim; direction++){

					if (userInputs.BC_list[starting_BC_list_index+component].var_BC_type[direction] == DIRICHLET){
						rigidBodyMode = false;
					}

				}
				// If the component has a rigid body mode, add it to the list
				if (rigidBodyMode == true){
					rigidBodyModeComponents.push_back(component);
				}
			}
		}
}

// Set constraints to pin the solution if there are no Dirichlet BCs for a component of a variable in an elliptic equation
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::setRigidBodyModeConstraints(const std::vector<int> rigidBodyModeComponents, ConstraintMatrix * constraints, const DoFHandler<dim>* dof_handler) const {

	if ( rigidBodyModeComponents.size() > 0 ){

		// Choose the point where the constraint will be placed. Must be the coordinates of a vertex.
		dealii::Point<dim> target_point; // default constructor places the point at the origin

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

#include "../../include/matrixFreePDE_template_instantiations.h"
