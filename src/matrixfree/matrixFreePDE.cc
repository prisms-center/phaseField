//constructor and destructor for matrixFreePDE class

#include "../../include/matrixFreePDE.h"


 //constructor
template <int dim, int degree>
 MatrixFreePDE<dim,degree>::MatrixFreePDE (userInputParameters<dim> _userInputs)
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
 currentOutput(0),
 computing_timer (pcout, TimerOutput::summary, TimerOutput::wall_times),
 first_integrated_var_output_complete(false)
 {
 }

 //destructor
 template <int dim, int degree>
 MatrixFreePDE<dim,degree>::~MatrixFreePDE ()
 {
   matrixFreeObject.clear();
   for(unsigned int iter=0; iter<fields.size(); iter++){
     //delete soltransSet[iter];
     delete locally_relevant_dofsSet[iter];
     delete constraintsDirichletSet[iter];
     delete dofHandlersSet[iter];
     delete FESet[iter];
     delete solutionSet[iter];
     delete residualSet[iter];
   }
 }


#include "../../include/matrixFreePDE_template_instantiations.h"
