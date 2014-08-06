//beta prime precipitate evolution implementation
//general headers
#include "../../include/dealIIheaders.h"

//precipitate problem headers
#include "parameters.h"
#include "../../src/CH.h"

//initial condition functions
//concentration initial conditions
template <int dim>
double InitialConditionC<dim>::value (const Point<dim> &p, const unsigned int /* component */) const
{
  //return the value of the initial concentration field at point p 
  return  0.6+ 0.2*(0.5 - (double)(std::rand() % 100 )/100.0);
}

//main
int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
      deallog.depth_console(0);
      CahnHilliardProblem<problemDIM> problem;
      problem.run ();
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
