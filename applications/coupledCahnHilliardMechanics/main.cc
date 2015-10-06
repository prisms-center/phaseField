//Coupled Cahn-Hilliard and Mechanics problem 
//general headers
#include "../../include/dealIIheaders.h"

//Coupled Cahn-Hilliard+Mechanics problem headers
#include "parameters.h"
#include "../../src/models/coupled/coupledCHMechanics.h"

//initial condition for concentration
template <int dim>
class InitialConditionC : public Function<dim>
{
public:
  InitialConditionC () : Function<dim>(1) {
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
void CoupledCahnHilliardMechanicsProblem<dim>::applyInitialConditions()
{
  unsigned int fieldIndex;
  //call initial condition function for c
  fieldIndex=this->getFieldIndex("c");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex],		\
			    InitialConditionC<dim>(),			\
			    *this->solutionSet[fieldIndex]);
  //set zero intial condition for mu
  fieldIndex=this->getFieldIndex("mu");
  *this->solutionSet[fieldIndex]=0.0;
  //set zero intial condition for u
  fieldIndex=this->getFieldIndex("u");
  *this->solutionSet[fieldIndex]=0.0;
}

//apply Dirchlet BC function
template <int dim>
void CoupledCahnHilliardMechanicsProblem<dim>::applyDirichletBCs(){
  //Set u=0 at all boundaries 
  VectorTools::interpolate_boundary_values (*this->dofHandlersSet[this->getFieldIndex("u")],\
					    0, ZeroFunction<dim>(dim), *(ConstraintMatrix*) \
					    this->constraintsSet[this->getFieldIndex("u")]);
}


//main
int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
      deallog.depth_console(0);
      CoupledCahnHilliardMechanicsProblem<problemDIM> problem;
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "mu"));
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "c"));
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
