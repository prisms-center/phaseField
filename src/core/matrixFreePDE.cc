// constructor and destructor for matrixFreePDE class

#include <core/matrixFreePDE.h>

// constructor
template <int dim, int degree>
MatrixFreePDE<dim, degree>::MatrixFreePDE(userInputParameters<dim> _userInputs)
  : Subscriptor()
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , userInputs(_userInputs)
  , var_attributes(_userInputs.var_attributes)
  , pp_attributes(_userInputs.pp_attributes)
  , triangulation(MPI_COMM_WORLD)
  , currentFieldIndex(0)
  , isTimeDependentBVP(false)
  , isEllipticBVP(false)
  , hasExplicitEquation(false)
  , hasNonExplicitEquation(false)
  , currentTime(0.0)
  , currentIncrement(0)
  , currentOutput(0)
  , currentCheckpoint(0)
  , current_grain_reassignment(0)
  , computing_timer(pcout, TimerOutput::summary, TimerOutput::wall_times)
  , first_integrated_var_output_complete(false)
  , AMR(_userInputs,
        triangulation,
        fields,
        solutionSet,
        soltransSet,
        FESet,
        dofHandlersSet_nonconst,
        constraintsDirichletSet,
        constraintsOtherSet)
{}

// destructor
template <int dim, int degree>
MatrixFreePDE<dim, degree>::~MatrixFreePDE()
{
  matrixFreeObject.clear();

  // Delete the pointers contained in several member variable vectors
  // The size of each of these must be checked individually in case an exception
  // is thrown as they are being initialized.
  for (const auto &locally_relevant_dofs : locally_relevant_dofsSet)
    {
      delete locally_relevant_dofs;
    }
  for (const auto &constraintsDirichlet : constraintsDirichletSet)
    {
      delete constraintsDirichlet;
    }
  for (const auto &soltrans : soltransSet)
    {
      delete soltrans;
    }
  for (const auto &dofHandlers : dofHandlersSet)
    {
      delete dofHandlers;
    }
  for (const auto &FE : FESet)
    {
      delete FE;
    }
  for (const auto &solution : solutionSet)
    {
      delete solution;
    }
  for (const auto &residual : residualSet)
    {
      delete residual;
    }
}

template class MatrixFreePDE<2, 1>;
template class MatrixFreePDE<3, 1>;

template class MatrixFreePDE<2, 2>;
template class MatrixFreePDE<3, 2>;

template class MatrixFreePDE<3, 3>;
template class MatrixFreePDE<2, 3>;

template class MatrixFreePDE<3, 4>;
template class MatrixFreePDE<2, 4>;

template class MatrixFreePDE<3, 5>;
template class MatrixFreePDE<2, 5>;

template class MatrixFreePDE<3, 6>;
template class MatrixFreePDE<2, 6>;