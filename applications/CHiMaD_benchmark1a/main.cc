// Cahn-Hilliard CHiMaD benchmark problem

// Header files
#include "parameters.h"
#include "../../include/dealIIheaders.h"
#include "../../include/typeDefs.h"
#include "../../include/model_variables.h"
#include "../../include/varBCs.h"
#include "../../include/loadInputs.h"
#include "../../include/initialConditions.h"
#include "../../include/matrixFreePDE.h"
#include "customPDE.h"
#include "equations.h"
#include "ICs_and_BCs.h"
#include "../../include/initialCondition_template_instantiations.h"
#include "../../src/loadInputs/loadInputs.cc" // Needs to be included because it contains needs access to the define macros in the preceding files

//main
int main (int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,dealii::numbers::invalid_unsigned_int);
  try
    {
	  dealii::deallog.depth_console(0);

	  userInputParameters userInputs;
	  userInputs.loadUserInput();

	  customPDE<problemDIM,finiteElementDegree> problem(userInputs);


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
