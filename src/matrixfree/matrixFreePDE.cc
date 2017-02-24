//constructor and destructor for matrixFreePDE class

#ifndef MATRIXFREEPDE_MATRIXFREE_H
#define MATRIXFREEPDE_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "init.cc"
#include "reinit.cc"
#include "initForTests.cc"
#include "refine.cc"
#include "invM.cc"
#include "computeLHS.cc"
#include "computeRHS.cc"
#include "modifyFields.cc"
#include "solve.cc"
#include "solveIncrement.cc"
#include "outputResults.cc"
#include "markBoundaries.cc"
#include "boundaryConditions.cc"
#include "initialConditions.cc"
#include "utilities.cc"
#include "calcFreeEnergy.cc"
#include "integrate_and_shift_field.cc"
#include "getOutputTimeSteps.cc"
#include "buildFields.cc"

 //constructor
 template <int dim>
 MatrixFreePDE<dim>::MatrixFreePDE (userInputParameters _userInputs)
 :
 Subscriptor(),
 userInputs(_userInputs),
 triangulation (MPI_COMM_WORLD),
 isTimeDependentBVP(false),
 isEllipticBVP(false),
 currentTime(0.0),
 currentIncrement(0),
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
