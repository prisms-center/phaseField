// init() method for MatrixFreePDE class
 
#ifndef INIT_MATRIXFREE_H
#define INIT_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

 //populate with fields and setup matrix free system
 template <int dim>
 void MatrixFreePDE<dim>::init(){
	 computing_timer.enter_section("matrixFreePDE: initialization");

	 //creating mesh

	 pcout << "creating problem mesh...\n";

     #if problemDIM==3
	 GridGenerator::subdivided_hyper_rectangle (triangulation, userInputs.subdivisions, Point<dim>(), Point<dim>(userInputs.domain_size[0],userInputs.domain_size[1],userInputs.domain_size[2]));
     #elif problemDIM==2
	 GridGenerator::subdivided_hyper_rectangle (triangulation, userInputs.subdivisions, Point<dim>(), Point<dim>(userInputs.domain_size[0],userInputs.domain_size[1]));
     #elif problemDIM==1
	 GridGenerator::subdivided_hyper_rectangle (triangulation, userInputs.subdivisions, Point<dim>(), Point<dim>(userInputs.domain_size[0]));
     #endif

	 // Mark boundaries for applying the boundary conditions
	 markBoundaries();

	 // Set which (if any) faces of the triangulation are periodic
	 setPeriodicity();

	 // Do the initial global refinement
	 triangulation.refine_global (userInputs.refine_factor);

	 // Write out the size of the computational domain and the total number of elements
	 pcout << "problem dimensions: " << userInputs.domain_size[0] << "x" << userInputs.domain_size[2] << "x" << userInputs.domain_size[3] << std::endl;
	 pcout << "number of elements: " << triangulation.n_global_active_cells() << std::endl;
	 pcout << std::endl;
  
	 // Setup system
	 pcout << "initializing matrix free object\n";
	 unsigned int totalDOFs=0;
	 for(typename std::vector<Field<dim> >::iterator it = fields.begin(); it != fields.end(); ++it){
		 currentFieldIndex=it->index;

		 char buffer[100];

		 //print to std::out
		 sprintf(buffer,"initializing finite element space P^%u for %9s:%6s field '%s'\n", \
				 userInputs.fe_degree,					\
			   (it->pdetype==PARABOLIC ? "PARABOLIC":"ELLIPTIC"),	\
			   (it->type==SCALAR ? "SCALAR":"VECTOR"),			\
			   it->name.c_str());
		 pcout << buffer;

		 // Check if any time dependent fields present (note: I should get rid of parabolicFieldIndex and ellipticFieldIndex, they only work if there is at max one of each)
		 if (it->pdetype==PARABOLIC){
			 isTimeDependentBVP=true;
			 parabolicFieldIndex=it->index;
		 }
		 else if (it->pdetype==ELLIPTIC){
			 isEllipticBVP=true;
			 ellipticFieldIndex=it->index;
		 }

		 //create FESystem
		 FESystem<dim>* fe;

		 if (it->type==SCALAR){
			 fe=new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(userInputs.fe_degree+1)),1);
		 }
		 else if (it->type==VECTOR){
			 fe=new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(userInputs.fe_degree+1)),dim);
		 }
		 else{
			 pcout << "\nmatrixFreePDE.h: unknown field type\n";
			 exit(-1);
		 }
		 FESet.push_back(fe);

		 //distribute DOFs
		 DoFHandler<dim>* dof_handler;

		 dof_handler=new DoFHandler<dim>(triangulation);
		 dofHandlersSet.push_back(dof_handler);
		 dofHandlersSet_nonconst.push_back(dof_handler);

		 dof_handler->distribute_dofs (*fe);
		 totalDOFs+=dof_handler->n_dofs();

		 // Extract locally_relevant_dofs
		 IndexSet* locally_relevant_dofs;

		 locally_relevant_dofs=new IndexSet;
		 locally_relevant_dofsSet.push_back(locally_relevant_dofs);
		 locally_relevant_dofsSet_nonconst.push_back(locally_relevant_dofs);

		 locally_relevant_dofs->clear();
		 DoFTools::extract_locally_relevant_dofs (*dof_handler, *locally_relevant_dofs);

		 // Create constraints
		 ConstraintMatrix *constraintsDirichlet, *constraintsOther;

		 constraintsDirichlet=new ConstraintMatrix; constraintsDirichletSet.push_back(constraintsDirichlet);
		 constraintsDirichletSet_nonconst.push_back(constraintsDirichlet);
		 constraintsOther=new ConstraintMatrix; constraintsOtherSet.push_back(constraintsOther);
		 constraintsOtherSet_nonconst.push_back(constraintsOther);
		 valuesDirichletSet.push_back(new std::map<dealii::types::global_dof_index, double>);

		 constraintsDirichlet->clear(); constraintsDirichlet->reinit(*locally_relevant_dofs);
		 constraintsOther->clear(); constraintsOther->reinit(*locally_relevant_dofs);

		 // Get hanging node constraints
		 DoFTools::make_hanging_node_constraints (*dof_handler, *constraintsOther);

		 // Add a constraint to fix the value at the origin to zero if all BCs are zero-derivative or periodic
		 std::vector<int> rigidBodyModeComponents;
		 getComponentsWithRigidBodyModes(rigidBodyModeComponents);
		 setRigidBodyModeConstraints(rigidBodyModeComponents,constraintsOther,dof_handler);

		 // Get constraints for periodic BCs
		 setPeriodicityConstraints(constraintsOther,dof_handler);

		 // Get constraints for Dirichlet BCs
		 applyDirichletBCs();

		 constraintsDirichlet->close();
		 constraintsOther->close();

		 // Store Dirichlet BC DOF's
		 valuesDirichletSet[it->index]->clear();
		 for (types::global_dof_index i=0; i<dof_handler->n_dofs(); i++){
			 if (locally_relevant_dofs->is_element(i)){
				 if (constraintsDirichlet->is_constrained(i)){
					 (*valuesDirichletSet[it->index])[i] = constraintsDirichlet->get_inhomogeneity(i);
				 }
			 }
		 }

		 sprintf(buffer, "field '%2s' DOF : %u (Constraint DOF : %u)\n", \
				 it->name.c_str(), dof_handler->n_dofs(), constraintsDirichlet->n_constraints());
		 pcout << buffer;
	 }
	 pcout << "total DOF : " << totalDOFs << std::endl;

	 // Setup the matrix free object
	 typename MatrixFree<dim,double>::AdditionalData additional_data;
	 additional_data.mpi_communicator = MPI_COMM_WORLD;
	 additional_data.tasks_parallel_scheme = MatrixFree<dim,double>::AdditionalData::partition_partition;
	 additional_data.mapping_update_flags = (update_values | update_gradients | update_JxW_values | update_quadrature_points);
	 QGaussLobatto<1> quadrature (userInputs.fe_degree+1);
	 matrixFreeObject.clear();
	 matrixFreeObject.reinit (dofHandlersSet, constraintsOtherSet, quadrature, additional_data);

	 bool dU_scalar_init = false;
	 bool dU_vector_init = false;
 
	 // Setup solution vectors
	 pcout << "initializing parallel::distributed residual and solution vectors\n";
	 for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
		 vectorType *U, *R;

		 U=new vectorType; R=new vectorType;
		 solutionSet.push_back(U); residualSet.push_back(R);
		 matrixFreeObject.initialize_dof_vector(*R,  fieldIndex); *R=0;

		 matrixFreeObject.initialize_dof_vector(*U,  fieldIndex); *U=0;

		 // Initializing temporary dU vector required for implicit solves of the elliptic equation.
		 // Assuming here that there is only one elliptic field in the problem (the main problem is if one is a scalar and the other is a vector, because then dU would need to be different sizes)
		 if (fields[fieldIndex].pdetype==ELLIPTIC){
			 if (fields[fieldIndex].type == SCALAR){
				 if (dU_scalar_init == false){
					 matrixFreeObject.initialize_dof_vector(dU_scalar,  fieldIndex);
					 dU_scalar_init = true;
				 }
			 }
			 else {
				 if (dU_vector_init == false){
					 matrixFreeObject.initialize_dof_vector(dU_vector,  fieldIndex);
					 dU_vector_init = true;
				 }
			 }
		 }
	 }
   
	 //check if time dependent BVP and compute invM
	 if (isTimeDependentBVP){
		 computeInvM();
	 }
   
	 // Apply the initial conditions to the solution vectors
	 // The initial conditions are re-applied below in the "adaptiveRefine" function so that the mesh can
	 // adapt based on the initial conditions.
	 applyInitialConditions();


	 // Create new solution transfer sets (needed for the "refineGrid" call, might be able to move this elsewhere)
	 soltransSet.clear();
	 for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
		 soltransSet.push_back(new parallel::distributed::SolutionTransfer<dim, vectorType>(*dofHandlersSet_nonconst[fieldIndex]));
	 }
   
	 // Ghost the solution vectors. Also apply the Dirichet BC's (if any) on the solution vectors
	 for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
		 constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
		 solutionSet[fieldIndex]->update_ghost_values();
	 }

	 // Check and perform adaptive mesh refinement, which reinitializes the system with the new mesh
	 adaptiveRefine(0);

	 computing_timer.exit_section("matrixFreePDE: initialization");
}

#endif 
