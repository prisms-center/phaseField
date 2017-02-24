//mesh refinement methods for MatrixFreePDE class
#ifndef REFINE_MATRIXFREE_H
#define REFINE_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//default implementation of adaptive mesh refinement 
template <int dim>
void MatrixFreePDE<dim>::adaptiveRefine(unsigned int currentIncrement){
#if hAdaptivity == true
	if ( (currentIncrement == 0) ){
		for (unsigned int remesh_index=0; remesh_index < (userInputs.max_refinement_level-userInputs.min_refinement_level); remesh_index++){
			this->reinit();
		}
	}
	else if ( (currentIncrement%skipRemeshingSteps==0) ){
			this->reinit();
		}
	#endif
}

//default implementation of adaptive mesh criterion
template <int dim>
void MatrixFreePDE<dim>::adaptiveRefineCriterion(){
  //Kelly error estimation criterion
  //estimate cell wise errors for mesh refinement
//#if hAdaptivity==true
//#ifdef adaptivityType
//#if adaptivityType=="KELLY"
//  Vector<float> estimated_error_per_cell (this->triangulation.n_locally_owned_active_cells());
//  KellyErrorEstimator<dim>::estimate (*this->dofHandlersSet_nonconst[refinementDOF],
//				      QGaussLobatto<dim-1>(finiteElementDegree+1),
//				      typename FunctionMap<dim>::type(),
//				      *this->solutionSet[refinementDOF],
//				      estimated_error_per_cell,
//				      ComponentMask(),
//				      0,
//				      1,
//				      this->triangulation.locally_owned_subdomain());
//  //flag cells for refinement
//  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction (this->triangulation,
//									    estimated_error_per_cell,
//									    topRefineFraction,
//									    bottomCoarsenFraction);
//#endif
//#endif
//#endif
#if hAdaptivity == true
	//Custom defined estimation criterion

	std::vector<std::vector<double> > errorOutV;


	QGauss<dim>  quadrature(userInputs.fe_degree+1);
	FEValues<dim> fe_values (*this->FESet[refine_criterion_fields[0]], quadrature, update_values);
	const unsigned int   num_quad_points = quadrature.size();

	std::vector<double> errorOut(num_quad_points);

	typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandlersSet_nonconst[refine_criterion_fields[0]]->begin_active(), endc = this->dofHandlersSet_nonconst[refine_criterion_fields[0]]->end();

	for (;cell!=endc; ++cell){
		if (cell->is_locally_owned()){
			fe_values.reinit (cell);

			for (unsigned int field_index=0; field_index<refine_criterion_fields.size(); field_index++){
				fe_values.get_function_values(*this->solutionSet[refine_criterion_fields[field_index]], errorOut);
				errorOutV.push_back(errorOut);
			}

			bool mark_refine = false;

			for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
				for (unsigned int field_index=0; field_index<refine_criterion_fields.size(); field_index++){
					if ((errorOutV[field_index][q_point]>refine_window_min[field_index]) && (errorOutV[field_index][q_point]<refine_window_max[field_index])){
						mark_refine = true;
						break;
					}
				}
			}

			errorOutV.clear();

//			fe_values.get_function_values(*this->solutionSet[refine_criterion_fields[0]], errorOut);
//
//			bool mark_refine = false;
//
//			for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
//				if ((errorOut[q_point]>refine_window_min[0]) && (errorOut[q_point]<refine_window_max[0])){
//					mark_refine = true;
//					break;
//				}
//			}


			if ( (mark_refine == true) ){
				cell->set_refine_flag();
			}
			else {
				cell->set_coarsen_flag();
			}
		}
	}
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
    for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(userInputs.max_refinement_level); cell != triangulation.end(); ++cell)
      cell->clear_refine_flag ();
  }

  //limit the minimal refinement depth of the mesh
  for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(userInputs.min_refinement_level); cell != triangulation.end(); ++cell){
      if (cell->level() <= userInputs.min_refinement_level ){
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
