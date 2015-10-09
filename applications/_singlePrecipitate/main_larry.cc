// Beta prime precipitate evolution implementation
// Code to calculate the steady-state morphology of a single precipitate
//general headers
#include "../../include/dealIIheaders.h"

//precipitate problem headers
#include "parameters.h"
#include "../../src/coupled_CH_AC_Mechanics.h"
 
// Set the initial condition for the concentration and the structural order parameter

//concentration initial conditions
template <int dim>
double InitialConditionC<dim>::value (const Point<dim> &p, const unsigned int /* component */) const
{
  //set result equal to the concentration initial condition 

  // Calculate the grid spacing
  double dz=spanZ/( (double)subdivisionsZ )/(std::pow(2.0,refineFactor));
  // Initialize the distance from the center of the domain to zero
  double r=0.0;
#if problemDIM==1
  r=p[0];
  return 0.005+0.5*(0.125-0.005)*(1-std::tanh((r-spanX/2.0)/(3*dx)));
#elif problemDIM==2
  //r=p.distance(Point<dim>(spanX/2.0,spanY/2.0));	
  r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/144.0+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0));
  return 0.5*(0.125-0.0)*(1.0-std::tanh((r-spanY/16.0)/(2.0*dy))) +0.03;
#elif problemDIM==3
  //r=p.distance(Point<dim>(spanX/2.0,spanY/2.0,spanZ/2.0));
  r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/4.0
  	+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/4.0
  	+(p.operator()(2)-spanZ/2.0)*(p.operator()(2)-spanZ/2.0)/1.0);
  return 0.5*(0.12-0.0)*(1.0-std::tanh((r-2.3811)/(1.0*dz))) +0.000;
#endif
}

//structural order parameter initial conditions
template <int dim>
double InitialConditionN<dim>::value (const Point<dim> &p, const unsigned int /* component */) const
{
  //set result equal to the structural order paramter initial condition
  double dz=spanZ/( (double)subdivisionsZ )/(std::pow(2.0,refineFactor));
  double r=0.0;
#if problemDIM==1
  r=p[0];
  return 0.5*(1.0-std::tanh((r-spanX/2.0)/(6.2*dx)));
#elif problemDIM==2
  if (index==0){
    //double r2=p.distance(Point<dim>(3*spanX/4.0,3*spanY/4.0));
    //r=std::min(r1,r2);
    //r=p.distance(Point<dim>(spanX/2.0,spanY/2.0));
	r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/144.0+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0));
	return 0.5*(1.0-std::tanh((r-spanY/16.0)/(2.0*dy)));

  }
  else if (index==1){
	return 0.0;  
//     double r1=p.distance(Point<dim>(3*spanX/4.0,spanY/4.0));
//     double r2=p.distance(Point<dim>(spanX/2.0,spanY/2.0));
//     r=std::min(r1,r2);
  }
  else if (index==2){
	return 0.0;  
//     r=p.distance(Point<dim>(spanX/4.0,3*spanY/4.0));
  }
  //return 0.5*(1.0-std::tanh((r-spanX/16.0)/(3*dx)));
#elif problemDIM==3
  if (index==0){
  //r=p.distance(Point<dim>(spanX/2.0,spanY/2.0,spanZ/2.0));
  //return 0.5*(1.0-std::tanh((r-spanX/8.0)/(3*dx)));
  r=sqrt((p.operator()(0)-spanX/2.0)*(p.operator()(0)-spanX/2.0)/4.0
  	+(p.operator()(1)-spanY/2.0)*(p.operator()(1)-spanY/2.0)/4.0
  	+(p.operator()(2)-spanZ/2.0)*(p.operator()(2)-spanZ/2.0)/1.0);
	return 0.5*(1.0-std::tanh((r-2.3811)/(1.0*dz)));
  }
  else if (index==1){
	return 0.0;  
  }
  else if (index==2){
	return 0.0;  
  }
  
#endif
}

//main
int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
      deallog.depth_console(0);
      //PrecipitateProblem<problemDIM> precipitateProblem;
      CoupledCHACMechanicsProblem<problemDIM> problem;
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
