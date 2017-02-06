// Coupled Cahn-Hilliard/Allen-Cahn example application

// Header files
#include "../../include/dealIIheaders.h"
#include "parameters.h"
#include "../../src/models/coupled/generalized_model.h"
#include <time.h>
#include <random>
#include "equations.h"
#include "ICs_and_BCs.h"
#include "../../src/models/coupled/generalized_model_functions.h"
#include "nucleation.h"


// =================================================================================
// MAIN
// =================================================================================
int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
	  deallog.depth_console(0);
	  generalizedProblem<problemDIM> problem;

      problem.setBCs();
      problem.buildFields();
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
