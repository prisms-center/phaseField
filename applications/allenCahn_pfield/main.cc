// Allen-Cahn example application

// Header files
#include "IntegrationTools/PField.hh"


#include "../../include/dealIIheaders.h"

#include "parameters.h"
#include "../../src/models/coupled/generalized_model.h"
#include "equations.h"
#include "ICs_and_BCs.h"
#include "../../src/models/coupled/generalized_model_functions.h"

//typedef PRISMS::PField<double*, double, 2> ScalarField2D;
//ScalarField2D conc;

//main
int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {

//	  typedef PRISMS::Body<double*, 2> Body2D;
//	  double coord[2];
//	  double tmp[2];
//	  Body2D body;
//	  body.read_vtk("solution-020000.0.vtu");
//	  //ScalarField2D &conc = body.find_scalar_field("n");
//	  conc = body.find_scalar_field("n");
//
//	  // Printing out a general point
//	  coord[0] = 5.128964;
//	  coord[1] = 6.12596472;
//	  std::cout << "Eta is: " << conc(coord) << std::endl;

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
