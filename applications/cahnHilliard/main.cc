//Cahn-Hilliard spinodal decomposition implementation
//general headers
#include "../../include/dealIIheaders.h"

//Cahn-Hilliard problem headers
#include "parameters.h"
#include "../../src/models/diffusion/CH.h"

//initial condition function for concentration
template <int dim>
class InitialConditionC : public Function<dim>
{
public:
  InitialConditionC () : Function<dim>(1) {
    //seeding the random number generator
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
void CahnHilliardProblem<dim>::applyInitialConditions()
{
  unsigned int fieldIndex;
  //call initial condition function for c
  fieldIndex=this->getFieldIndex("c");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex],	\
			    InitialConditionC<dim>(),			\
			    *this->solutionSet[fieldIndex]);
  //set initial condition for mu
  fieldIndex=this->getFieldIndex("mu");
  *(this->residualSet[fieldIndex])=0.0;
  this->matrixFreeObject.cell_loop (&CahnHilliardProblem<dim>::getRHS, this, this->residualSet, this->solutionSet);
  //sove for mu from initial condition for c
  for (unsigned int dof=0; dof<this->solutionSet[fieldIndex]->local_size(); ++dof){
    this->solutionSet[fieldIndex]->local_element(dof)= this->invM.local_element(dof)*this->residualSet[fieldIndex]->local_element(dof);
  }
}

//main
int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
      deallog.depth_console(0);
      CahnHilliardProblem<problemDIM> problem;
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "mu"));
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "c"));
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
