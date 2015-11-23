//Coupled Cahn-Hilliard and Allen-Cahn implementation with nucleation
//general headers
#include "../../include/dealIIheaders.h"
//coupled Cahn-Hilliard and Allen-Cahn problem headers
#include "parameters.h"
#include "../../src/models/coupled/coupledCHACMechanics.h"
#include <time.h>

//initial condition function for concentration
template <int dim>
class InitialConditionC : public Function<dim>
{
public:
  InitialConditionC () : Function<dim>(1) {
    //std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
    std::srand(time(NULL));
  }
  double value (const Point<dim> &p, const unsigned int component = 0) const
  {
    //return the value of the initial concentration field at point p 
    return avg_Nd;
  }
};

//initial condition function for the order parameters
template <int dim>
class InitialConditionN : public Function<dim>
{
public:
  unsigned int index;
  InitialConditionN (const unsigned int _index) : Function<dim>(1), index(_index) {
    std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const Point<dim> &p, const unsigned int component = 0) const
  {
    //set result equal to the structural order paramter initial condition
	  double dx=spanX/( (double)subdivisionsX )/std::pow(2.0,refineFactor);
	  double dy=spanY/( (double)subdivisionsY )/std::pow(2.0,refineFactor);
	  double dz=spanZ/( (double)subdivisionsZ )/std::pow(2.0,refineFactor);
    double r=0.0;
	#if problemDIM==1
    r=p.operator()(0);
	  return 0.5*(1.0-std::tanh((r-spanX/16.0)/(0.1*dx)));
	#elif problemDIM==2
	  if (index==1){
		  //r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0));
		  //return 0.5*(1.0-std::tanh((r-spanX/20.0)/(2.0*dx)));
		  return 0.0;
	  }
	  else if (index==2){
		return 0.0;
	  }
	  else if (index==3){
		return 0.0;
	  }
	  //return 0.5*(1.0-std::tanh((r-spanX/16.0)/(3*dx)));
	#elif problemDIM==3
	  if (index==1){
	  //r=p.distance(Point<dim>(spanX/2.0,spanY/2.0,spanZ/2.0));
	  //return 0.5*(1.0-std::tanh((r-spanX/8.0)/(3*dx)));
		  r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)
		  		+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)
		  		+(p.operator()(2)-spanZ/2.0)*(p.operator()(2)-spanZ/2.0));
		return 0.5*(1.0-std::tanh((r-2.3811)/(1.0*dz)));
	  }
	  else if (index==2){
		return 0.0;
	  }
	  else if (index==3){
		return 0.0;
	  }
	#endif
	return 0.0;
  }
};

//apply initial conditions
template <int dim>
void CoupledCHACMechanicsProblem<dim>::applyInitialConditions()
{
  unsigned int fieldIndex;
  
  //call initial condition function for c
  fieldIndex=this->getFieldIndex("c");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex],		\
			    InitialConditionC<dim>(),			\
			    *this->solutionSet[fieldIndex]);
  //call initial condition function for structural order parameters
  fieldIndex=this->getFieldIndex("n1");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex], InitialConditionN<dim>(1), *this->solutionSet[fieldIndex]);
  fieldIndex=this->getFieldIndex("n2");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex], InitialConditionN<dim>(2), *this->solutionSet[fieldIndex]);
  fieldIndex=this->getFieldIndex("n3");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex], InitialConditionN<dim>(3), *this->solutionSet[fieldIndex]);
  //set zero intial condition for u
  fieldIndex=this->getFieldIndex("u");
  *this->solutionSet[fieldIndex]=0.0;
}

//apply Dirchlet BC function
template <int dim>
void CoupledCHACMechanicsProblem<dim>::applyDirichletBCs(){
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
      CoupledCHACMechanicsProblem<problemDIM> problem;
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "c"));
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "n1"));
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "n2"));
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "n3"));
      problem.fields.push_back(Field<problemDIM>(VECTOR,  ELLIPTIC, "u"));
      problem.init(); 
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
