//methods to apply initial conditions 

#ifndef INITIALCONDITIONS_MATRIXFREE_H
#define INITIALCONDITIONS_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//methods to apply initial conditions
template <int dim>
void MatrixFreePDE<dim>::applyInitialConditions(){
  pcout << "applying the default zero initial condition on all fields\n";
  //default method to apply zero initial conditions on all fields
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    VectorTools::interpolate (*this->dofHandlersSet[fieldIndex],	\
			      ZeroFunction<dim>(fields[fieldIndex].numComponents), \
			      *this->solutionSet[fieldIndex]);
  }
}
#endif
