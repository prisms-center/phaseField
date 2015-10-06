//Fickian diffusion problem
//general headers
#include "../../include/dealIIheaders.h"

//Fickian diffusion problem headers
#include "parameters.h"
#include "../../src/models/diffusion/Fickian.h"

//apply initial conditions
template <int dim>
void FickianProblem<dim>::applyInitialConditions()
{
  unsigned int fieldIndex;
  fieldIndex=this->getFieldIndex("c");
  //set c=0
  *this->solutionSet[fieldIndex]=0.0;
}

//main
int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
      deallog.depth_console(0);
      FickianProblem<problemDIM> problem;
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "c"));
      problem.init (); 
      problem.solve();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  
  return 0;
}
