//Coupled Cahn-Hilliard and Allen-Cahn implementation
//general headers
#include "../../../include/dealIIheaders.h"

//coupled Cahn-Hilliard and Allen-Cahn problem headers
#include "parameters.h"
#include "../../../src/models/diffusion/coupledCHAC.h"
 
//initial condition function for concentration
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
    double dx=spanX/((double) subdivisionsX)/std::pow(2.0,refineFactor);
    double r=0.0;
#if problemDIM==1
    r=p[0];
    return 0.005+0.5*(0.125-0.005)*(1-std::tanh((r-spanX/2.0)/(3*dx)));
#elif problemDIM==2
    r=p.distance(Point<dim>(spanX/2.0,spanY/2.0));
    return 0.005+0.5*(0.125-0.005)*(1-std::tanh((r-spanX/8.0)/(3*dx)));
#elif problemDIM==3
    r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/x_denom
         		+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/y_denom
         		+(p.operator()(2)-spanZ/2.0)*(p.operator()(2)-spanZ/2.0)/z_denom);
    return 0.5*(0.12-avg_c)*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff))) + avg_c;

#endif
  }
};

//initial condition function for order parameter
template <int dim>
class InitialConditionN : public Function<dim>
{
public:
  InitialConditionN () : Function<dim>(1) {
    std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const Point<dim> &p, const unsigned int component = 0) const
  {
    //return the value of the initial order parameter field at point p 
	  double dx=spanX/((double) subdivisionsX)/std::pow(2.0,refineFactor);
    double r=0.0;
#if problemDIM==1
  r=p[0];
  return 0.5*(1.0-std::tanh((r-spanX/2.0)/(6.2*dx)));
#elif problemDIM==2
  r=p.distance(Point<dim>(spanX/2.0,spanY/2.0));
  return 0.5*(1.0-std::tanh((r-spanX/8.0)/(3*dx)));
#elif problemDIM==3

  r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/x_denom
     		+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/y_denom
     		+(p.operator()(2)-spanZ/2.0)*(p.operator()(2)-spanZ/2.0)/z_denom);
     return 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));

#endif
  }
};

//apply initial conditions
template <int dim>
void CoupledCHACProblem<dim>::applyInitialConditions()
{
  unsigned int fieldIndex;
  //call initial condition function for c
  fieldIndex=this->getFieldIndex("c");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex],		\
			    InitialConditionC<dim>(),			\
			    *this->solutionSet[fieldIndex]);
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
      CoupledCHACProblem<problemDIM> problem;
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "n"));
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
