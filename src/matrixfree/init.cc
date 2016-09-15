//init() method for MatrixFreePDE class
 
#ifndef INIT_MATRIXFREE_H
#define INIT_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

 //populate with fields and setup matrix free system
 template <int dim>
 void MatrixFreePDE<dim>::init(unsigned int iter){
   //void MatrixFreePDE<dim>::init(std::vector<Field<dim> >& _fields){
   computing_timer.enter_section("matrixFreePDE: initialization"); 

   if (iter==0){
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

     triangulation.refine_global (refineFactor);
     //write out extends
     pcout << "problem dimensions: " << spanX << "x" << spanY << "x" << spanZ << std::endl;
     pcout << "number of elements: " << triangulation.n_global_active_cells() << std::endl;
     pcout << std::endl;
  
     //mark boundaries for applying Dirichlet boundary conditons
     markBoundaries();
   }
   else{
     refineGrid();
   }
     
   //setup system
   pcout << "initializing matrix free object\n";
   unsigned int totalDOFs=0;
   for(typename std::vector<Field<dim> >::iterator it = fields.begin(); it != fields.end(); ++it){
     char buffer[100];
     if (iter==0){
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
     }

     
     //create FESystem
     FESystem<dim>* fe;
     //
     if (iter==0){
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
     }
     else{
       fe=FESet.at(it->index);
     }
     
     //distribute DOFs
     DoFHandler<dim>* dof_handler;
     if (iter==0){
       dof_handler=new DoFHandler<dim>(triangulation);
       dofHandlersSet.push_back(dof_handler);
       dofHandlersSet2.push_back(dof_handler); 
     }
     else{
       dof_handler=dofHandlersSet2.at(it->index);
     }
     dof_handler->distribute_dofs (*fe);
     totalDOFs+=dof_handler->n_dofs();

     //extract locally_relevant_dofs
     IndexSet* locally_relevant_dofs;
     if (iter==0){
       locally_relevant_dofs=new IndexSet;
       locally_relevant_dofsSet.push_back(locally_relevant_dofs);
       locally_relevant_dofsSet2.push_back(locally_relevant_dofs);
     }
     else{
       locally_relevant_dofs=locally_relevant_dofsSet2.at(it->index);
     }
     locally_relevant_dofs->clear();
     DoFTools::extract_locally_relevant_dofs (*dof_handler, *locally_relevant_dofs);

     //create constraints
     ConstraintMatrix *constraints, *constraintsHangingNodes;
     if (iter==0){
       constraints=new ConstraintMatrix; constraintsSet.push_back(constraints);
       constraintsSet2.push_back(constraints);
       constraintsHangingNodes=new ConstraintMatrix; constraintsHangingNodesSet.push_back(constraintsHangingNodes);
       constraintsHangingNodesSet2.push_back(constraintsHangingNodes);
       valuesDirichletSet.push_back(new std::map<dealii::types::global_dof_index, double>);
     }
     else{
       constraints=constraintsSet2.at(it->index);
       constraintsHangingNodes=constraintsHangingNodesSet2.at(it->index);
     }
     constraints->clear(); constraints->reinit(*locally_relevant_dofs);
     constraintsHangingNodes->clear(); constraintsHangingNodes->reinit(*locally_relevant_dofs);
     DoFTools::make_hanging_node_constraints (*dof_handler, *constraintsHangingNodes);
     
     //apply Dirichlet BC's
     currentFieldIndex=it->index;
     applyDirichletBCs();

     constraints->close();
     constraintsHangingNodes->close();

     //store Dirichlet BC DOF's
     valuesDirichletSet[it->index]->clear();
     for (types::global_dof_index i=0; i<dof_handler->n_dofs(); i++){
       if (locally_relevant_dofs->is_element(i)){
	 if (constraints->is_constrained(i)){
	   (*valuesDirichletSet[it->index])[i] = constraints->get_inhomogeneity(i);
	 }
       }
     }

     // Periodic BCs
     std::vector<dealii::GridTools::PeriodicFacePair<typename dealii::parallel::distributed::Triangulation<dim>::cell_iterator> > periodicity_vector;
    // dealii::GridTools::collect_periodic_faces(triangulation, 0, 1, 0,
      //                                 periodicity_vector, Tensor<1, dim>());


//     //periodic BCs from CHiMaD
//
//     //make periodic BC's
//        std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> > periodicity_vector;
//        for (int i=0; i<dim; ++i){
//          GridTools::collect_periodic_faces(triangulation, /*b_id1*/ 2*i, /*b_id2*/ 2*i+1,
//     				       /*direction*/ i, periodicity_vector);
//        }
//        triangulation.add_periodicity(periodicity_vector);
//        std::cout << "periodic facepairs: " << periodicity_vector.size() << std::endl;
//

//	  std::vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator> > periodicity_vector1;
//	  for (int i=0; i<dim; ++i){
//		GridTools::collect_periodic_faces(*dof_handler, /*b_id1*/ 2*i, /*b_id2*/ 2*i+1,
//					 /*direction*/ i, periodicity_vector1);
//	  }
//	  DoFTools::make_periodicity_constraints<DoFHandler<dim> >(periodicity_vector1, *constraints);


     sprintf(buffer, "field '%2s' DOF : %u (Constraint DOF : %u)\n", \
	     it->name.c_str(), dof_handler->n_dofs(), constraints->n_constraints());
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
   matrixFreeObject.reinit (dofHandlersSet, constraintsHangingNodesSet, quadrature, additional_data);
 
   //setup problem vectors
   pcout << "initializing parallel::distributed residual and solution vectors\n";
   for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
     vectorType *U, *R;
     if (iter==0){
       U=new vectorType; R=new vectorType;
       solutionSet.push_back(U); residualSet.push_back(R); 
       matrixFreeObject.initialize_dof_vector(*R,  fieldIndex); *R=0;
     }
     else{
       U=solutionSet.at(fieldIndex); 
     }
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
   
   //apply initial conditions if iter=0, else transfer solution from previous refined mesh
   if (iter==0){
     applyInitialConditions();
   }
   else{
     for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
       //interpolate and clear used solution transfer sets
       soltransSet[fieldIndex]->interpolate(*solutionSet[fieldIndex]);
       delete soltransSet[fieldIndex];
       
       //reset residual vector
       vectorType *R=residualSet.at(fieldIndex);
       matrixFreeObject.initialize_dof_vector(*R,  fieldIndex); *R=0;
     }
   }

   //apply Dirichlet BC's
   if (isEllipticBVP){
     constraintsSet[ellipticFieldIndex]->distribute(*solutionSet.at(ellipticFieldIndex));
     constraintsHangingNodesSet[ellipticFieldIndex]->distribute(*solutionSet.at(ellipticFieldIndex));
   }
   
   //create new solution transfer sets
   soltransSet.clear();
   for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
     soltransSet.push_back(new parallel::distributed::SolutionTransfer<dim, vectorType>(*dofHandlersSet2[fieldIndex]));
   }
   
   //Ghost the solution vectors. Also apply the Dirichet BC's (if any) on the solution vectors 
   for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
     constraintsSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
     solutionSet[fieldIndex]->update_ghost_values();
   } 

   computing_timer.exit_section("matrixFreePDE: initialization");  
}

#endif 
