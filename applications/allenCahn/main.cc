//Allen-Cahn order parameter evolution implementation
//general headers
//general headers
#include "../../include/dealIIheaders.h"

//Allen-Hilliard problem headers
#include "parameters.h"
#include "residuals.h"
#include "../../src/models/diffusion/AC.h"
#include "ICs_and_BCs.h"

//main
int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
      deallog.depth_console(0);
      AllenCahnProblem<problemDIM> problem;
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "n"));
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

