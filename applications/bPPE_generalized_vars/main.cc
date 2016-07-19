//Coupled Cahn-Hilliard, Allen-Cahn and Mechanics problem
//general headers
#include "../../include/dealIIheaders.h"

//Coupled Cahn-Hilliard+Allen-Cahn+Mechanics problem headers
#include "parameters.h"
#include "../../src/models/coupled/coupledCHACMechanics_generalized.h"
#include "ICs_and_BCs.h"
#include "../../src/models/coupled/generalized_model_functions.h"
#include "residuals.h"


//main
int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
	  deallog.depth_console(0);
      CoupledCHACMechanicsProblem<problemDIM> problem;

      problem.setBCs();

      // Build each of the fields in the system
      for (unsigned int i=0; i<num_var; i++){
    	  if (var_type[i] == "SCALAR"){
    		  if (var_eq_type[i] == "ELLIPTIC"){
    			  problem.fields.push_back(Field<problemDIM>(SCALAR, ELLIPTIC, var_name[i]));
    		  }
    		  else if (var_eq_type[i] == "PARABOLIC"){
    			  problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, var_name[i]));
    		  }
    		  else{
    			  // Need to change to throw an exception
    			  std::cerr << "Error: Equation type must be ELLIPTIC or PARABOLIC " << std::endl;
    		  }
    	  }
    	  else if (var_type[i] == "VECTOR"){
    		  if (var_eq_type[i] == "ELLIPTIC"){
    			  problem.fields.push_back(Field<problemDIM>(VECTOR, ELLIPTIC, var_name[i]));
    		  }
    		  else if (var_eq_type[i] == "PARABOLIC"){
    			  problem.fields.push_back(Field<problemDIM>(VECTOR, PARABOLIC, var_name[i]));
    		  }
    		  else{
    			  // Need to change to throw an exception
    			  std::cerr << "Error: Variable type must be SCALAR or VECTOR " << std::endl;
    		  }
    	  }
      }

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
