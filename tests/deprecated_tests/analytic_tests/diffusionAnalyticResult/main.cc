// Coupled Cahn-Hilliard, Allen-Cahn and Mechanics problem
// general headers
#include "../../../include/dealIIheaders.h"

// Coupled Cahn-Hilliard+Allen-Cahn+Mechanics problem headers
#include "ICs_and_BCs.h"
#include "parameters.h"
#include "residuals.h"

#include "../../../src/models/coupled/generalized_model.h"
#include "../../../src/models/coupled/generalized_model_functions.h"

// main
int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                      argv,
                                                      numbers::invalid_unsigned_int);
  try
    {
      deallog.depth_console(0);
      generalizedProblem<problemDIM> problem;

      problem.setBCs();
      problem.buildFields();
      problem.init();
      problem.solve();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------" << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------" << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------" << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------" << std::endl;
      return 1;
    }

  return 0;
}
