//init() method for MatrixFreePDE class

#ifndef INIT_MATRIXFREE_H
#define INIT_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

 //populate with fields and setup matrix free system
 template <int dim>
 void MatrixFreePDE<dim>::init(){
 //void MatrixFreePDE<dim>::init(std::vector<Field<dim> >& _fields){
   computing_timer.enter_section("matrixFreePDE: initialization"); 

   //creating mesh
   pcout << "creating problem mesh\n";
 #if problemDIM==3
   GridGenerator::hyper_rectangle (triangulation, Point<dim>(), Point<dim>(spanX,spanY,spanZ));
#elif problemDIM==2
   GridGenerator::hyper_rectangle (triangulation, Point<dim>(), Point<dim>(spanX,spanY));
 #elif problemDIM==1
   GridGenerator::hyper_rectangle (triangulation, Point<dim>(), Point<dim>(spanX));
 #endif
   triangulation.refine_global (refineFactor);
   pcout << "number of elements: " << triangulation.n_global_active_cells() << std::endl;

   //setup system
   unsigned int totalDOFs=0;
   for(typename std::vector<Field<dim> >::iterator it = fields.begin(); it != fields.end(); ++it){
     //print to std::out
     char buffer[100];
     sprintf(buffer,"initializing finite element space P^%u for %6s:%9s field '%s'\n", \
	     finiteElementDegree,					\
	     (it->type==SCALAR ? "SCALAR":"VECTOR"),			\
	     (it->pdetype==PARABOLIC ? "PARABOLIC":"ELLIPTIC"), it->name.c_str());
     pcout << buffer;

     //check if any time dependent fields present
     if (it->pdetype==PARABOLIC){
       isTimeDependentBVP=true;
     }

     //create FESystem
     FESystem<dim>* fe;
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
     DoFHandler<dim>* dof_handler=new DoFHandler<dim>(triangulation);
     dof_handler->distribute_dofs (*fe);
     dofHandlersSet.push_back(dof_handler); 
     totalDOFs+=dof_handler->n_dofs();

     //extract locally_relevant_dofs
     IndexSet* locally_relevant_dofs=new IndexSet;
     DoFTools::extract_locally_relevant_dofs (*dof_handler, *locally_relevant_dofs);
     locally_relevant_dofsSet.push_back(locally_relevant_dofs);

     //create constraints
     ConstraintMatrix* constraints=new ConstraintMatrix;
     constraints->clear();
     constraints->reinit(*locally_relevant_dofs);
     DoFTools::make_hanging_node_constraints (*dof_handler, *constraints);
     //apply zero Dirichlet BC's for ELLIPTIC fields. This is just the
     //default and can be changed later in the specific BVP
     //implementation
     if (it->pdetype==ELLIPTIC){
       VectorTools::interpolate_boundary_values (*dof_handler, 0, ZeroFunction<dim>(it->numComponents), *constraints);
     }
     constraints->close();  
     constraintsSet.push_back(constraints);
   }
   pcout << "number of degrees of freedom: " << totalDOFs << std::endl;

   //setup the matrix free object
   typename MatrixFree<dim,double>::AdditionalData additional_data;
   additional_data.mpi_communicator = MPI_COMM_WORLD;
   additional_data.tasks_parallel_scheme = MatrixFree<dim,double>::AdditionalData::partition_partition;
   additional_data.mapping_update_flags = (update_values | update_gradients | update_JxW_values);
   QGaussLobatto<1> quadrature (finiteElementDegree+1);
   matrixFreeObject.reinit (dofHandlersSet, constraintsSet, quadrature, additional_data);
   pcout << "completed initialization of the matrix free object\n";
 
   //setup problem vectors
   pcout << "initializing parallel::distributed residual and solution vectors\n";
   for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
     vectorType* U=new vectorType;
     vectorType* R=new vectorType;
     matrixFreeObject.initialize_dof_vector(*U,  fieldIndex);
     matrixFreeObject.initialize_dof_vector(*R,  fieldIndex);
     *U=0; solutionSet.push_back(U);
     *R=0; residualSet.push_back(R);
   }

   //check if time dependent BVP and compute invM
   if (isTimeDependentBVP){
     computeInvM();
   }
   computing_timer.exit_section("matrixFreePDE: initialization");  
}

#endif
