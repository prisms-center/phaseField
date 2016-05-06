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
  //set initial condition for mu
  fieldIndex=this->getFieldIndex("mu");
  *(this->residualSet[fieldIndex])=0.0;
  this->matrixFreeObject.cell_loop (&CoupledCahnHilliardMechanicsProblem<dim>::getRHS, this, this->residualSet, this->solutionSet);
  //sove for mu from initial condition for c
  for (unsigned int dof=0; dof<this->solutionSet[fieldIndex]->local_size(); ++dof){
    this->solutionSet[fieldIndex]->local_element(dof)= this->invM.local_element(dof)*this->residualSet[fieldIndex]->local_element(dof);
  }
  //set zero initial condition for u
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

