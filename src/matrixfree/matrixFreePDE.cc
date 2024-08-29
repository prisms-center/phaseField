// constructor and destructor for matrixFreePDE class

#include "../../include/matrixFreePDE.h"

// constructor
template <int dim, int degree>
MatrixFreePDE<dim, degree>::MatrixFreePDE(userInputParameters<dim> _userInputs)
  : Subscriptor()
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , userInputs(_userInputs)
  , triangulation(MPI_COMM_WORLD)
  , currentFieldIndex(0)
  , isTimeDependentBVP(false)
  , isEllipticBVP(false)
  , hasExplicitEquation(false)
  , hasNonExplicitEquation(false)
  , parabolicFieldIndex(0)
  , ellipticFieldIndex(0)
  , currentTime(0.0)
  , currentIncrement(0)
  , currentOutput(0)
  , currentCheckpoint(0)
  , current_grain_reassignment(0)
  , computing_timer(pcout, TimerOutput::summary, TimerOutput::wall_times)
  , first_integrated_var_output_complete(false)
{}

// destructor
template <int dim, int degree>
MatrixFreePDE<dim, degree>::~MatrixFreePDE()
{
  matrixFreeObject.clear();

  // Delete the pointers contained in several member variable vectors
  // The size of each of these must be checked individually in case an exception
  // is thrown as they are being initialized.
  for (unsigned int iter = 0; iter < locally_relevant_dofsSet.size(); iter++)
    {
      delete locally_relevant_dofsSet[iter];
    }
  for (unsigned int iter = 0; iter < constraintsDirichletSet.size(); iter++)
    {
      delete constraintsDirichletSet[iter];
    }
  for (unsigned int iter = 0; iter < soltransSet.size(); iter++)
    {
      delete soltransSet[iter];
    }
  for (unsigned int iter = 0; iter < dofHandlersSet.size(); iter++)
    {
      delete dofHandlersSet[iter];
    }
  for (unsigned int iter = 0; iter < FESet.size(); iter++)
    {
      delete FESet[iter];
    }
  for (unsigned int iter = 0; iter < solutionSet.size(); iter++)
    {
      delete solutionSet[iter];
    }
  for (unsigned int iter = 0; iter < residualSet.size(); iter++)
    {
      delete residualSet[iter];
    }
}

#include "../../include/matrixFreePDE_template_instantiations.h"
