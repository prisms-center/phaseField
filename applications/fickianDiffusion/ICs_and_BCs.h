//apply initial conditions
template <int dim>
void FickianProblem<dim>::applyInitialConditions()
{
  unsigned int fieldIndex;
  fieldIndex=this->getFieldIndex("c");
  //set c=0
  *this->solutionSet[fieldIndex]=0.0;
}
