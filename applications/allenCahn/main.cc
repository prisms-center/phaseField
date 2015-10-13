//Allen-Cahn order parameter evolution implementation
//general headers
//general headers
#include "../../include/dealIIheaders.h"

//Allen-Hilliard problem headers
#include "parameters.h"
#include "../../src/models/diffusion/AC.h"

//initial condition function for the order parameter
template <int dim>
class InitialConditionN : public Function<dim>
{
public:
  InitialConditionN () : Function<dim>(1) {
    std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const Point<dim> &p, const unsigned int component = 0) const
  {
    //return the value of the initial concentration field at point p 
    return  0.5+ 0.2*(0.5 - (double)(std::rand() % 100 )/100.0);
  }
};

//apply initial conditions
template <int dim>
void AllenCahnProblem<dim>::applyInitialConditions()
{
  unsigned int fieldIndex;
  //call initial condition function for n
  fieldIndex=this->getFieldIndex("n");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex],		\
			    InitialConditionN<dim>(),			\
			    *this->solutionSet[fieldIndex]);
}

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

