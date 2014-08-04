//beta prime precipitate evolution implementation
//general headers
#include "../../include/dealIIheaders.h"

//precipitate problem headers
#include "parameters.h"
#include "../../src/coupled_CH_AC.h"

//initial condition functions
//concentration initial conditions
template <int dim>
double InitialConditionC<dim>::value (const Point<dim> &p, const unsigned int /* component */) const
{
  //set result equal to the concentration initial condition 
  double result =  0.02 + 1.0e-3*(2*(0.5 - (double)(std::rand() % 100 )/100.0));
  return result;
}

//structural order parameter initial conditions
template <int dim>
double InitialConditionN<dim>::value (const Point<dim> &p, const unsigned int /* component */) const
{
  //set result equal to the structural order paramter initial condition
  double result = 0.0;
  if (p.distance(Point<dim>(50.0,50.0))<16.0) result=1.0;
  return result;
}

//main
int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
      deallog.depth_console(0);
      PrecipitateProblem<problemDIM> precipitateProblem;
      precipitateProblem.run ();
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
