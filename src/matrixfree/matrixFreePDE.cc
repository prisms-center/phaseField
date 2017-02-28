//constructor and destructor for matrixFreePDE class

#ifndef MATRIXFREEPDE_MATRIXFREE_H
#define MATRIXFREEPDE_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../include/matrixFreePDE.h"


 //constructor
template <int dim, int degree>
 MatrixFreePDE<dim,degree>::MatrixFreePDE (userInputParameters _userInputs)
 :
 Subscriptor(),
 pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0),
 userInputs(_userInputs),
 triangulation (MPI_COMM_WORLD),
 currentFieldIndex(0),
 isTimeDependentBVP(false),
 isEllipticBVP(false),
 parabolicFieldIndex(0),
 ellipticFieldIndex(0),
 currentTime(0.0),
 currentIncrement(0),
 computing_timer (pcout, TimerOutput::summary, TimerOutput::wall_times),
 energy(0.0)
 {
 }

 //destructor
 template <int dim, int degree>
 MatrixFreePDE<dim,degree>::~MatrixFreePDE ()
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


#ifndef MATRIXFREEPDE_TEMPLATE_INSTANTIATION
#define MATRIXFREEPDE_TEMPLATE_INSTANTIATION
template class MatrixFreePDE<2,1>;
template class MatrixFreePDE<3,1>;
#endif



#endif
