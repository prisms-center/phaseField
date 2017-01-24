////init() method for MatrixFreePDE class
//
//#ifndef INIT_MATRIXFREE_H
//#define INIT_MATRIXFREE_H
////this source file is temporarily treated as a header file (hence
////#ifndef's) till library packaging scheme is finalized
//
// //populate with fields and setup matrix free system
// template <int dim>
// void MatrixFreePDE<dim>::init(unsigned int iter){
//   //void MatrixFreePDE<dim>::init(std::vector<Field<dim> >& _fields){
//   computing_timer.enter_section("matrixFreePDE: initialization");
//
//   if (iter==0){
//     //creating mesh
//     std::vector<unsigned int> subdivisions;
//     subdivisions.push_back(subdivisionsX);
//     if (dim>1){
//       subdivisions.push_back(subdivisionsY);
//       if (dim>2){
//	 subdivisions.push_back(subdivisionsZ);
//       }
//     }
//
//     pcout << "creating problem mesh...\n";
//
//#ifndef t_domain
//#define t_domain false
//#endif
//#if t_domain == true
//     std::vector< unsigned int > repetitions(2);
//        repetitions[0]=15;//5*refineFactor;
//        repetitions[1]=3;//1*refineFactor;
//        Triangulation<2> tria1, tria2;
//        GridGenerator::subdivided_hyper_rectangle (tria1, repetitions, Point<dim>(-40,100), Point<dim>(60,120));
//        repetitions[0]=3;//*refineFactor;
//        repetitions[1]=15;//*refineFactor;
//        GridGenerator::subdivided_hyper_rectangle (tria2, repetitions, Point<dim>(0,0), Point<dim>(20,100));
//        GridGenerator::merge_triangulations (tria1, tria2, triangulation);
//#else
//	#if problemDIM==3
//     GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions, Point<dim>(), Point<dim>(spanX,spanY,spanZ));
//	#elif problemDIM==2
//     GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions, Point<dim>(), Point<dim>(spanX,spanY));
//	#elif problemDIM==1
//     GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions, Point<dim>(), Point<dim>(spanX));
//	#endif
//
//
//     // Mark boundaries for applying the boundary conditions
//     markBoundaries();
//
//     // Set which (if any) faces of the triangulation are periodic
//    setPeriodicity();
//
//#endif
//     // Do the initial global refinement
//     triangulation.refine_global (refineFactor);
//
//
//
//     // Write out the size of the computational domain and the total number of elements
//     pcout << "problem dimensions: " << spanX << "x" << spanY << "x" << spanZ << std::endl;
//     pcout << "number of elements: " << triangulation.n_global_active_cells() << std::endl;
//     pcout << std::endl;
//
//
//   }
//   else{
//     refineGrid();
//   }
//
//   //setup system
//   pcout << "initializing matrix free object\n";
//   unsigned int totalDOFs=0;
//   for(typename std::vector<Field<dim> >::iterator it = fields.begin(); it != fields.end(); ++it){
//	   currentFieldIndex=it->index;
//
//	   char buffer[100];
//	   if (iter==0){
//		   //print to std::out
//		   sprintf(buffer,"initializing finite element space P^%u for %9s:%6s field '%s'\n", \
//				   finiteElementDegree,					\
//				   (it->pdetype==PARABOLIC ? "PARABOLIC":"ELLIPTIC"),	\
//				   (it->type==SCALAR ? "SCALAR":"VECTOR"),			\
//				   it->name.c_str());
//		   pcout << buffer;
//		   //check if any time dependent fields present
//		   if (it->pdetype==PARABOLIC){
//			   isTimeDependentBVP=true;
//			   parabolicFieldIndex=it->index;
//		   }
//		   else if (it->pdetype==ELLIPTIC){
//			   isEllipticBVP=true;
//			   ellipticFieldIndex=it->index;
//		   }
//	   }
//
//
//	   //create FESystem
//	   FESystem<dim>* fe;
//	   //
//	   if (iter==0){
//		   if (it->type==SCALAR){
//			   fe=new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(finiteElementDegree+1)),1);
//		   }
//		   else if (it->type==VECTOR){
//			   fe=new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(finiteElementDegree+1)),dim);
//		   }
//		   else{
//			   pcout << "\nmatrixFreePDE.h: unknown field type\n";
//			   exit(-1);
//		   }
//		   FESet.push_back(fe);
//	   }
//	   else{
//		   fe=FESet.at(it->index);
//	   }
//
//	   //distribute DOFs
//	   DoFHandler<dim>* dof_handler;
//	   if (iter==0){
//		   dof_handler=new DoFHandler<dim>(triangulation);
//		   dofHandlersSet.push_back(dof_handler);
//		   dofHandlersSet_nonconst.push_back(dof_handler);
//	   }
//	   else{
//		   dof_handler=dofHandlersSet_nonconst.at(it->index);
//	   }
//	   dof_handler->distribute_dofs (*fe);
//	   totalDOFs+=dof_handler->n_dofs();
//
//	   //extract locally_relevant_dofs
//	   IndexSet* locally_relevant_dofs;
//	   if (iter==0){
//		   locally_relevant_dofs=new IndexSet;
//		   locally_relevant_dofsSet.push_back(locally_relevant_dofs);
//		   locally_relevant_dofsSet_nonconst.push_back(locally_relevant_dofs);
//	   }
//	   else{
//		   locally_relevant_dofs=locally_relevant_dofsSet_nonconst.at(it->index);
//	   }
//	   locally_relevant_dofs->clear();
//	   DoFTools::extract_locally_relevant_dofs (*dof_handler, *locally_relevant_dofs);
//
//	   //create constraints
//	   ConstraintMatrix *constraintsDirichlet, *constraintsOther;
//
//	   if (iter==0){
//		   constraintsDirichlet=new ConstraintMatrix; constraintsDirichletSet.push_back(constraintsDirichlet);
//		   constraintsDirichletSet_nonconst.push_back(constraintsDirichlet);
//		   constraintsOther=new ConstraintMatrix; constraintsOtherSet.push_back(constraintsOther);
//		   constraintsOtherSet_nonconst.push_back(constraintsOther);
//		   valuesDirichletSet.push_back(new std::map<dealii::types::global_dof_index, double>);
//	   }
//	   else{
//		   constraintsDirichlet=constraintsDirichletSet_nonconst.at(it->index);
//		   constraintsOther=constraintsOtherSet_nonconst.at(it->index);
//	   }
//	   constraintsDirichlet->clear(); constraintsDirichlet->reinit(*locally_relevant_dofs);
//	   constraintsOther->clear(); constraintsOther->reinit(*locally_relevant_dofs);
//	   DoFTools::make_hanging_node_constraints (*dof_handler, *constraintsOther);
//
//	   // Add a constraint to fix the value at the origin to zero if all BCs are zero-derivative or periodic
//	   std::vector<int> rigidBodyModeComponents;
//	   getComponentsWithRigidBodyModes(rigidBodyModeComponents);
//	   setRigidBodyModeConstraints(rigidBodyModeComponents,constraintsOther,dof_handler);
//
//	   // Apply periodic BCs
//	   setPeriodicityConstraints(constraintsOther,dof_handler);
//
//
//	   // Apply Dirichlet BCs
//	   applyDirichletBCs();
//
//	   constraintsDirichlet->close();
//	   constraintsOther->close();
//
//	   //store Dirichlet BC DOF's
//	   valuesDirichletSet[it->index]->clear();
//	   for (types::global_dof_index i=0; i<dof_handler->n_dofs(); i++){
//		   if (locally_relevant_dofs->is_element(i)){
//			   if (constraintsDirichlet->is_constrained(i)){
//				   (*valuesDirichletSet[it->index])[i] = constraintsDirichlet->get_inhomogeneity(i);
//			   }
//		   }
//	   }
//
//	   sprintf(buffer, "field '%2s' DOF : %u (Constraint DOF : %u)\n", \
//			   it->name.c_str(), dof_handler->n_dofs(), constraintsDirichlet->n_constraints());
//	   pcout << buffer;
//   }
//   pcout << "total DOF : " << totalDOFs << std::endl;
//
//   //setup the matrix free object
//   typename MatrixFree<dim,double>::AdditionalData additional_data;
//   additional_data.mpi_communicator = MPI_COMM_WORLD;
//   additional_data.tasks_parallel_scheme = MatrixFree<dim,double>::AdditionalData::partition_partition;
//   additional_data.mapping_update_flags = (update_values | update_gradients | update_JxW_values | update_quadrature_points);
//   QGaussLobatto<1> quadrature (finiteElementDegree+1);
//   num_quadrature_points=std::pow(quadrature.size(),dim);
//   matrixFreeObject.clear();
//   matrixFreeObject.reinit (dofHandlersSet, constraintsOtherSet, quadrature, additional_data);
//
//   //setup problem vectors
//   pcout << "initializing parallel::distributed residual and solution vectors\n";
//   for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
//     vectorType *U, *R;
//     if (iter==0){
//       U=new vectorType; R=new vectorType;
//       solutionSet.push_back(U); residualSet.push_back(R);
//       matrixFreeObject.initialize_dof_vector(*R,  fieldIndex); *R=0;
//     }
//     else{
//       U=solutionSet.at(fieldIndex);
//     }
//     matrixFreeObject.initialize_dof_vector(*U,  fieldIndex); *U=0;
//
//     //initializing temporary dU vector required for implicit solves of the elliptic equation.
//     //Assuming here that there is only one elliptic field in the problem
//     if (fields[fieldIndex].pdetype==ELLIPTIC){
//    	 matrixFreeObject.initialize_dof_vector(dU,  fieldIndex);
//     }
//   }
//
//   //check if time dependent BVP and compute invM
//   if (isTimeDependentBVP){
//     computeInvM();
//   }
//
//   //apply initial conditions if iter=0, else transfer solution from previous refined mesh
//   if (iter==0){
//	   applyInitialConditions();
//   }
//   else{
//     for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
//       //interpolate and clear used solution transfer sets
//       soltransSet[fieldIndex]->interpolate(*solutionSet[fieldIndex]);
//       delete soltransSet[fieldIndex];
//
//       //reset residual vector
//       vectorType *R=residualSet.at(fieldIndex);
//       matrixFreeObject.initialize_dof_vector(*R,  fieldIndex); *R=0;
//     }
//   }
//
//   //apply Dirichlet BC's
//   if (isEllipticBVP){
//     constraintsDirichletSet[ellipticFieldIndex]->distribute(*solutionSet.at(ellipticFieldIndex));
//     constraintsOtherSet[ellipticFieldIndex]->distribute(*solutionSet.at(ellipticFieldIndex));
//   }
//
//   //create new solution transfer sets
//   soltransSet.clear();
//   for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
//     soltransSet.push_back(new parallel::distributed::SolutionTransfer<dim, vectorType>(*dofHandlersSet_nonconst[fieldIndex]));
//   }
//
//   //apply initial conditions if iter=0, else transfer solution from previous refined mesh
//   if (iter==0){
//	   computing_timer.exit_section("matrixFreePDE: initialization");
//	   adaptiveRefine(0);
//	   computing_timer.enter_section("matrixFreePDE: initialization");
//   	   applyInitialConditions();
//   }
//
//   //Ghost the solution vectors. Also apply the Dirichet BC's (if any) on the solution vectors
//   for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
//     constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
//     solutionSet[fieldIndex]->update_ghost_values();
//   }
//
//   computing_timer.exit_section("matrixFreePDE: initialization");
//}
//
//#endif

//reinit() method for MatrixFreePDE class
 
#ifndef INIT_MATRIXFREE_H
#define INIT_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

 //populate with fields and setup matrix free system
 template <int dim>
 void MatrixFreePDE<dim>::init(unsigned int iter){
   //void MatrixFreePDE<dim>::init(std::vector<Field<dim> >& _fields){
   computing_timer.enter_section("matrixFreePDE: initialization"); 

   //creating mesh
   std::vector<unsigned int> subdivisions;
   subdivisions.push_back(subdivisionsX);
   if (dim>1){
	   subdivisions.push_back(subdivisionsY);
	   if (dim>2){
		   subdivisions.push_back(subdivisionsZ);
	   }
   }

   pcout << "creating problem mesh...\n";

   #if problemDIM==3
   GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions, Point<dim>(), Point<dim>(spanX,spanY,spanZ));
   #elif problemDIM==2
   GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions, Point<dim>(), Point<dim>(spanX,spanY));
   #elif problemDIM==1
   GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions, Point<dim>(), Point<dim>(spanX));
   #endif

   // Mark boundaries for applying the boundary conditions
   markBoundaries();

   // Set which (if any) faces of the triangulation are periodic
   setPeriodicity();

   // Do the initial global refinement
   triangulation.refine_global (refineFactor);

   // Write out the size of the computational domain and the total number of elements
   pcout << "problem dimensions: " << spanX << "x" << spanY << "x" << spanZ << std::endl;
   pcout << "number of elements: " << triangulation.n_global_active_cells() << std::endl;
   pcout << std::endl;
  
   //setup system
   pcout << "initializing matrix free object\n";
   unsigned int totalDOFs=0;
   for(typename std::vector<Field<dim> >::iterator it = fields.begin(); it != fields.end(); ++it){
	   currentFieldIndex=it->index;

	   char buffer[100];

	   //print to std::out
	   sprintf(buffer,"initializing finite element space P^%u for %9s:%6s field '%s'\n", \
			   finiteElementDegree,					\
			   (it->pdetype==PARABOLIC ? "PARABOLIC":"ELLIPTIC"),	\
			   (it->type==SCALAR ? "SCALAR":"VECTOR"),			\
			   it->name.c_str());
	   pcout << buffer;
	   //check if any time dependent fields present
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
	   //

	   if (it->type==SCALAR){
		   fe=new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(finiteElementDegree+1)),1);
	   }
	   else if (it->type==VECTOR){
		   fe=new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(finiteElementDegree+1)),dim);
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

	   //extract locally_relevant_dofs
	   IndexSet* locally_relevant_dofs;

	   locally_relevant_dofs=new IndexSet;
	   locally_relevant_dofsSet.push_back(locally_relevant_dofs);
	   locally_relevant_dofsSet_nonconst.push_back(locally_relevant_dofs);

	   locally_relevant_dofs->clear();
	   DoFTools::extract_locally_relevant_dofs (*dof_handler, *locally_relevant_dofs);

	   //create constraints
	   ConstraintMatrix *constraintsDirichlet, *constraintsOther;

	   constraintsDirichlet=new ConstraintMatrix; constraintsDirichletSet.push_back(constraintsDirichlet);
	   constraintsDirichletSet_nonconst.push_back(constraintsDirichlet);
	   constraintsOther=new ConstraintMatrix; constraintsOtherSet.push_back(constraintsOther);
	   constraintsOtherSet_nonconst.push_back(constraintsOther);
	   valuesDirichletSet.push_back(new std::map<dealii::types::global_dof_index, double>);

	   constraintsDirichlet->clear(); constraintsDirichlet->reinit(*locally_relevant_dofs);
	   constraintsOther->clear(); constraintsOther->reinit(*locally_relevant_dofs);
	   DoFTools::make_hanging_node_constraints (*dof_handler, *constraintsOther);

	   // Add a constraint to fix the value at the origin to zero if all BCs are zero-derivative or periodic
	   std::vector<int> rigidBodyModeComponents;
	   getComponentsWithRigidBodyModes(rigidBodyModeComponents);
	   setRigidBodyModeConstraints(rigidBodyModeComponents,constraintsOther,dof_handler);

	   // Apply periodic BCs
	   setPeriodicityConstraints(constraintsOther,dof_handler);


	   // Apply Dirichlet BCs
	   applyDirichletBCs();

	   constraintsDirichlet->close();
	   constraintsOther->close();

	   //store Dirichlet BC DOF's
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

   //setup the matrix free object
   typename MatrixFree<dim,double>::AdditionalData additional_data;
   additional_data.mpi_communicator = MPI_COMM_WORLD;
   additional_data.tasks_parallel_scheme = MatrixFree<dim,double>::AdditionalData::partition_partition;
   additional_data.mapping_update_flags = (update_values | update_gradients | update_JxW_values | update_quadrature_points);
   QGaussLobatto<1> quadrature (finiteElementDegree+1);
   num_quadrature_points=std::pow(quadrature.size(),dim);
   matrixFreeObject.clear();
   matrixFreeObject.reinit (dofHandlersSet, constraintsOtherSet, quadrature, additional_data);
 
   //setup problem vectors
   pcout << "initializing parallel::distributed residual and solution vectors\n";
   for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
     vectorType *U, *R;

       U=new vectorType; R=new vectorType;
       solutionSet.push_back(U); residualSet.push_back(R); 
       matrixFreeObject.initialize_dof_vector(*R,  fieldIndex); *R=0;

     matrixFreeObject.initialize_dof_vector(*U,  fieldIndex); *U=0;
     
     //initializing temporary dU vector required for implicit solves of the elliptic equation.
     //Assuming here that there is only one elliptic field in the problem
     if (fields[fieldIndex].pdetype==ELLIPTIC){
    	 matrixFreeObject.initialize_dof_vector(dU,  fieldIndex);
     }
   }
   
   //check if time dependent BVP and compute invM
   if (isTimeDependentBVP){
     computeInvM();
   }
   
   // Apply the initial conditions.
   // The initial conditions are re-applied below in the "adaptiveRefine" function so that the mesh can
   // adapt based on the initial conditions.
   applyInitialConditions();

   // Apply Dirichlet BC's
//   if (isEllipticBVP){
//     constraintsDirichletSet[ellipticFieldIndex]->distribute(*solutionSet.at(ellipticFieldIndex));
//     constraintsOtherSet[ellipticFieldIndex]->distribute(*solutionSet.at(ellipticFieldIndex));
//   }
   
   //create new solution transfer sets (Code seg faults without this, why is it needed? Should only be needed if transferring between meshes)
   // Needed for the "refineGrid" call
   soltransSet.clear();
   for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
     soltransSet.push_back(new parallel::distributed::SolutionTransfer<dim, vectorType>(*dofHandlersSet_nonconst[fieldIndex]));
   }
   
   //Ghost the solution vectors. Also apply the Dirichet BC's (if any) on the solution vectors 
   for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
     constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
     solutionSet[fieldIndex]->update_ghost_values();
   }

   //check and perform adaptive mesh refinement
   computing_timer.exit_section("matrixFreePDE: initialization");
   adaptiveRefine(0);
   computing_timer.enter_section("matrixFreePDE: initialization");
   std::cout << "after refine (in init)!" << std::endl;

   // If adaptivity is turned off, apply initial conditions, Dirichlet constraints and ghost the solution vectors.
   // If adaptivity is turned on, these operations were performed in the "reinit" in the "adaptiveRefine" call.
//   if (hAdaptivity==false){
//
//	   applyInitialConditions();
//
//	   // Ghost the solution vectors. Also apply the Dirichet BC's (if any) on the solution vectors
//	   for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
//		   constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
//		   solutionSet[fieldIndex]->update_ghost_values();
//	   }
//   }

   computing_timer.exit_section("matrixFreePDE: initialization");  
}

#endif 
