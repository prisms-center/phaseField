//initForTests() method for MatrixFreePDE class

#ifndef INITFORTESTS_MATRIXFREE_H
#define INITFORTESTS_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

 template <int dim>
 void MatrixFreePDE<dim>::initForTests(){
   //creating mesh
   std::vector<unsigned int> subdivisions;
   subdivisions.push_back(10);
   if (dim>1){
     subdivisions.push_back(10);
     if (dim>2){
       subdivisions.push_back(10);
     }
   }
     
#if problemDIM==3
   GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions, Point<dim>(), Point<dim>(1,1,1));
#elif problemDIM==2
   GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions, Point<dim>(), Point<dim>(1,1));
#elif problemDIM==1
   GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions, Point<dim>(), Point<dim>(1));
#endif

   //setup system
   //create FESystem
   FESystem<dim>* fe;
   fe=new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(finiteElementDegree+1)),1);
   FESet.push_back(fe);

   //distribute DOFs
   DoFHandler<dim>* dof_handler;
   dof_handler=new DoFHandler<dim>(triangulation);
   dofHandlersSet.push_back(dof_handler);
   dofHandlersSet_nonconst.push_back(dof_handler);
   dof_handler->distribute_dofs (*fe);
   
   //extract locally_relevant_dofs
   IndexSet* locally_relevant_dofs;
   locally_relevant_dofs=new IndexSet;
   locally_relevant_dofsSet.push_back(locally_relevant_dofs);
   locally_relevant_dofsSet_nonconst.push_back(locally_relevant_dofs);
   locally_relevant_dofs->clear();
   DoFTools::extract_locally_relevant_dofs (*dof_handler, *locally_relevant_dofs);

   //create constraints
   ConstraintMatrix *constraintsOther;
   constraintsOther=new ConstraintMatrix; constraintsOtherSet.push_back(constraintsOther);
   constraintsOtherSet_nonconst.push_back(constraintsOther);
   constraintsOther->clear(); constraintsOther->reinit(*locally_relevant_dofs);
   DoFTools::make_hanging_node_constraints (*dof_handler, *constraintsOther);
   constraintsOther->close();

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
   vectorType *U, *R;
   U=new vectorType; R=new vectorType;
   solutionSet.push_back(U); residualSet.push_back(R); 
   matrixFreeObject.initialize_dof_vector(*R,  0); *R=0;
   matrixFreeObject.initialize_dof_vector(*U,  0); *U=0;

}

#endif 
