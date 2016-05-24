//mesh refinement methods for MatrixFreePDE class

#ifndef REFINE_MATRIXFREE_H
#define REFINE_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

template <int dim>
void MatrixFreePDE<dim>::refineGrid (){
  //estimate cell wise errors for mesh refinement
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  KellyErrorEstimator<dim>::estimate (*dofHandlersSet2[refinementDOF],
				      QGauss<dim-1>(finiteElementDegree+1),
				      typename FunctionMap<dim>::type(),
				      *solutionSet[refinementDOF],
				      estimated_error_per_cell);
  //flag cells for refinement
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
									    estimated_error_per_cell,
									    topRefineFraction,
									    bottomCoarsenFraction);
  //limit the maximal refinement depth of the mesh
  pcout << "levels: " << triangulation.n_levels() << "/n";
  if (triangulation.n_levels() > maxRefinementLevel){
    for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(maxRefinementLevel); cell != triangulation.end(); ++cell)
      cell->clear_refine_flag ();
  }

  //prepare and refine
  triangulation.prepare_coarsening_and_refinement();
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    (*residualSet[fieldIndex])=(*solutionSet[fieldIndex]);
    soltransSet[fieldIndex]->prepare_for_coarsening_and_refinement(*residualSet[fieldIndex]);
  }
  triangulation.execute_coarsening_and_refinement();
}


//initialize adaptive mesh refinement
template <int dim>
void MatrixFreePDE<dim>::refineMesh(unsigned int _currentIncrement){
  init(_currentIncrement-1);
}

//default implementation of adaptive mesh refinement 
template <int dim>
void MatrixFreePDE<dim>::adaptiveRefine(unsigned int currentIncrement){
}



#endif 
