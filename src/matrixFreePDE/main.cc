//Cahn-Hilliard spinodal decomposition implementation
//general headers
#include "../../include/dealIIheaders.h"

//Cahn-Hilliard problem headers
#include "parameters.h"
#include "matrixFreePDE.h"


//main
int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
      deallog.depth_console(0);
      MatrixFreePDE<problemDIM> problem;

      //
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "c"));
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "n"));
      problem.fields.push_back(Field<problemDIM>(VECTOR, ELLIPTIC, "u"));
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "v"));
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
