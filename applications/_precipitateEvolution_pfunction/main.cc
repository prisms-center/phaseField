// Beta prime precipitate evolution implementation
// Code to calculate the steady-state morphology of a single precipitate
//general headers
#include "../../include/dealIIheaders.h"

//Coupled Cahn-Hilliard+Allen-Cahn+Mechanics problem headers

// PLibrary includes
#include "parameters.h"
#include "../../src/models/coupled/coupledCHACMechanics.h"
#include "../_precipitateEvolution_pfunction/ICs_and_BCs.h"
#include "../_precipitateEvolution_pfunction/parameters.h"
#include "../_precipitateEvolution_pfunction/PLibrary/PLibrary.cc"
#include "../_precipitateEvolution_pfunction/PLibrary/PLibrary.hh"
#include "../_precipitateEvolution_pfunction/residuals.h"


//main
int main (int argc, char **argv)
{

	// Load variables from the PLibrary
	std::string homo_free_energy_funcname = "pfunct_faV";
	std::string mobility_funcname = "pfunct_McV";
	PRISMS::PLibrary::checkout(homo_free_energy_funcname, pfunct_faV);
	PRISMS::PLibrary::checkout("pfunct_fbV", pfunct_fbV);
	PRISMS::PLibrary::checkout(mobility_funcname, pfunct_McV);
	PRISMS::PLibrary::checkout("pfunct_Mn1V", pfunct_Mn1V);
	PRISMS::PLibrary::checkout("pfunct_Mn2V", pfunct_Mn2V);
	PRISMS::PLibrary::checkout("pfunct_Mn3V", pfunct_Mn3V);

	Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
	  deallog.depth_console(0);
	  CoupledCHACMechanicsProblem<problemDIM> problem;
	  problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "c"));
	  problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "n1"));
	  problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "n2"));
	  problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "n3"));
	  problem.fields.push_back(Field<problemDIM>(VECTOR,  ELLIPTIC, "u"));
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


