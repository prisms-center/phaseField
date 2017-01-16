//constructor and destructor for matrixFreePDE class

#ifndef MATRIXFREEPDE_MATRIXFREE_H
#define MATRIXFREEPDE_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

 //constructor
 template <int dim>
 MatrixFreePDE<dim>::MatrixFreePDE ()
 :
 Subscriptor(),
 triangulation (MPI_COMM_WORLD),
 isTimeDependentBVP(false),
 isEllipticBVP(false),
 dtValue(0.0),
 currentTime(0.0),
 finalTime(0.0),
 currentIncrement(0),
 totalIncrements(1),
 pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0),
 computing_timer (pcout, TimerOutput::summary, TimerOutput::wall_times)
 {

 }

 //destructor
 template <int dim>
 MatrixFreePDE<dim>::~MatrixFreePDE ()
 {
   matrixFreeObject.clear();
   for(unsigned int iter=0; iter<fields.size(); iter++){
     delete soltransSet[iter];
     delete locally_relevant_dofsSet[iter];
     delete constraintsDirichletSet[iter];
     delete dofHandlersSet[iter];
     delete FESet[iter];
     delete solutionSet[iter];
     delete residualSet[iter];
   } 
 }

#endif
