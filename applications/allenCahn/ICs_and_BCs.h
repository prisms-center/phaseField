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

