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
    r=p.distance(Point<dim>(spanX/2.0,spanY/2.0,spanZ/2.0));
    return 0.005+0.5*(0.125-0.005)*(1-std::tanh((r-spanX/8.0)/(3*dx)));
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
  r=p.distance(Point<dim>(spanX/2.0,spanY/2.0,spanZ/2.0));
  return 0.5*(1.0-std::tanh((r-spanX/8.0)/(3*dx)));
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
