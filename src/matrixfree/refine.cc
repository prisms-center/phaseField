//mesh refinement methods for MatrixFreePDE class
#ifndef REFINE_MATRIXFREE_H
#define REFINE_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//default implementation of adaptive mesh refinement 
template <int dim>
void MatrixFreePDE<dim>::adaptiveRefine(unsigned int currentIncrement){
}

//default implementation of adaptive mesh criterion
template <int dim>
void MatrixFreePDE<dim>::adaptiveRefineCriterion(){
  //Kelly error estimation criterion
  //estimate cell wise errors for mesh refinement
#if hAdaptivity==true 
#ifdef adaptivityType
#if adaptivityType=="KELLY"
  Vector<float> estimated_error_per_cell (this->triangulation.n_locally_owned_active_cells());
  KellyErrorEstimator<dim>::estimate (*this->dofHandlersSet_nonconst[refinementDOF],
				      QGaussLobatto<dim-1>(finiteElementDegree+1),
				      typename FunctionMap<dim>::type(),
				      *this->solutionSet[refinementDOF],
				      estimated_error_per_cell,
				      ComponentMask(),
				      0,
				      1,
				      this->triangulation.locally_owned_subdomain());
  //flag cells for refinement
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction (this->triangulation,
									    estimated_error_per_cell,
									    topRefineFraction,
									    bottomCoarsenFraction);
#endif
#endif
#endif
}


//refine grid method
template <int dim>
void MatrixFreePDE<dim>::refineGrid (){
#if hAdaptivity==true 
  //call refinement criterion for adaptivity
  adaptiveRefineCriterion();
  
  //limit the maximal refinement depth of the mesh
  pcout << "Current mesh refinement level: " << triangulation.n_levels() << "\n";
  if ( triangulation.n_levels() > maxRefinementLevel ){
    for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(maxRefinementLevel); cell != triangulation.end(); ++cell)
      cell->clear_refine_flag ();
  }

  //limit the minimal refinement depth of the mesh
  for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(minRefinementLevel); cell != triangulation.end(); ++cell){
      if (cell->level() <= minRefinementLevel ){
    	  cell->clear_coarsen_flag ();
      }
  }


  //prepare and refine
  triangulation.prepare_coarsening_and_refinement();
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    (*residualSet[fieldIndex])=(*solutionSet[fieldIndex]);
    soltransSet[fieldIndex]->prepare_for_coarsening_and_refinement(*residualSet[fieldIndex]);
  }
  triangulation.execute_coarsening_and_refinement();
#endif
}


#endif 
