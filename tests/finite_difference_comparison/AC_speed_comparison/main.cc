//Allen-Cahn order parameter evolution implementation
//general headers
//general headers
#include "../../../include/dealIIheaders.h"

//Allen-Hilliard problem headers
#include "parameters.h"
#include "../../../src/models/diffusion/AC.h"

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
	  //set result equal to the structural order paramter initial condition
	 	  double dx=spanX/( (double)subdivisionsX )/std::pow(2.0,refineFactor);
	 	  double dy=spanY/( (double)subdivisionsY )/std::pow(2.0,refineFactor);
	 	  double dz=spanZ/( (double)subdivisionsZ )/std::pow(2.0,refineFactor);
	     double r=0.0;
	 	#if problemDIM==1
	     r=p.operator()(0);
	 	  return 0.5*(1.0-std::tanh((r-spanX/16.0)/(0.1*dx)));
	 	#elif problemDIM==2

	 		  r=sqrt((p.operator()(0))*(p.operator()(0))/x_denom
	 		  		  		+(p.operator()(1))*(p.operator()(1))/y_denom);
	 		  return 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));
	 	  //return 0.5*(1.0-std::tanh((r-spanX/16.0)/(3*dx)));
	 	#elif problemDIM==3

	 	  //r=p.distance(Point<dim>(spanX/2.0,spanY/2.0,spanZ/2.0));
	 	  //return 0.5*(1.0-std::tanh((r-spanX/8.0)/(3*dx)));

	 		r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)
	 		  		+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)*4.0
	 		  		+(p.operator()(2)-spanZ/2.0)*(p.operator()(2)-spanZ/2.0));
	 		return 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));

	 		// planar interface
	 		//r=sqrt((p.operator()(2)-spanZ/2.0)*(p.operator()(2)-spanZ/2.0));
	 		//return 0.5*(1.0-std::tanh((r)/(initial_interface_coeff)));
	 		//return 0.5*(1.0-std::tanh((r-initial_radius)/(initial_interface_coeff)));
	 	#endif
	 	return 0.0;
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

