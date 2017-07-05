//mesh refinement methods for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

//default implementation of adaptive mesh refinement
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::adaptiveRefine(unsigned int currentIncrement){
if (userInputs.h_adaptivity == true){
	if ( (currentIncrement == 0) ){
		computing_timer.enter_section("matrixFreePDE: AMR");
		unsigned int numDoF_preremesh = totalDOFs;
		for (unsigned int remesh_index=0; remesh_index < (userInputs.max_refinement_level-userInputs.min_refinement_level); remesh_index++){

			adaptiveRefineCriterion();
			refineGrid();
			reinit();

			// If the mesh hasn't changed from the previous cycle, stop remeshing
			if (totalDOFs == numDoF_preremesh) break;
			numDoF_preremesh = totalDOFs;
		}
		computing_timer.exit_section("matrixFreePDE: AMR");
	}
	else if ( (currentIncrement%userInputs.skip_remeshing_steps==0) ){

		computing_timer.enter_section("matrixFreePDE: AMR");

		// Apply constraints before remeshing
		for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
			constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
			constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
			solutionSet[fieldIndex]->update_ghost_values();
		}
		adaptiveRefineCriterion();
		refineGrid();
		reinit();
		computing_timer.exit_section("matrixFreePDE: AMR");
	}
}
}

//default implementation of adaptive mesh criterion
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::adaptiveRefineCriterion(){
  //Kelly error estimation criterion
  //estimate cell wise errors for mesh refinement
//#if hAdaptivity==true
//#ifdef adaptivityType
//#if adaptivityType=="KELLY"
//  Vector<float> estimated_error_per_cell (triangulation.n_locally_owned_active_cells());
//  KellyErrorEstimator<dim>::estimate (*dofHandlersSet_nonconst[refinementDOF],
//				      QGaussLobatto<dim-1>(degree+1),
//				      typename FunctionMap<dim>::type(),
//				      *solutionSet[refinementDOF],
//				      estimated_error_per_cell,
//				      ComponentMask(),
//				      0,
//				      1,
//				      triangulation.locally_owned_subdomain());
//  //flag cells for refinement
//  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
//									    estimated_error_per_cell,
//									    topRefineFraction,
//									    bottomCoarsenFraction);
//#endif
//#endif
//#endif

//Custom defined estimation criterion

std::vector<std::vector<double> > errorOutV;


QGaussLobatto<dim>  quadrature(degree+1);
FEValues<dim> fe_values (*FESet[userInputs.refine_criterion_fields[0]], quadrature, update_values);
const unsigned int num_quad_points = quadrature.size();

std::vector<double> errorOut(num_quad_points);

typename DoFHandler<dim>::active_cell_iterator cell = dofHandlersSet_nonconst[userInputs.refine_criterion_fields[0]]->begin_active(), endc = dofHandlersSet_nonconst[userInputs.refine_criterion_fields[0]]->end();

typename parallel::distributed::Triangulation<dim>::active_cell_iterator t_cell = triangulation.begin_active();

for (;cell!=endc; ++cell){
	if (cell->is_locally_owned()){
		fe_values.reinit (cell);

		for (unsigned int field_index=0; field_index<userInputs.refine_criterion_fields.size(); field_index++){
			fe_values.get_function_values(*solutionSet[userInputs.refine_criterion_fields[field_index]], errorOut);
			errorOutV.push_back(errorOut);
		}

		bool mark_refine = false;

		for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
			for (unsigned int field_index=0; field_index<userInputs.refine_criterion_fields.size(); field_index++){
				if ((errorOutV[field_index][q_point]>userInputs.refine_window_min[field_index]) && (errorOutV[field_index][q_point]<userInputs.refine_window_max[field_index])){
					mark_refine = true;
					break;
				}
			}
		}

		errorOutV.clear();

		//limit the maximal and minimal refinement depth of the mesh
		unsigned int current_level = t_cell->level();

		if ( (mark_refine && current_level < userInputs.max_refinement_level) ){
			cell->set_refine_flag();
		}
		else if (!mark_refine && current_level > userInputs.min_refinement_level) {
			cell->set_coarsen_flag();
		}

	}
	++t_cell;
}

}


//refine grid method
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::refineGrid (){

//prepare and refine
triangulation.prepare_coarsening_and_refinement();
for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
	// The following lines were from an earlier version.
	// residualSet is cleared in reinit(), so I don't see the reason for the pointer assignment
	// Changing to the new version has no impact on the solution.
	//(*residualSet[fieldIndex])=(*solutionSet[fieldIndex]);
	//soltransSet[fieldIndex]->prepare_for_coarsening_and_refinement(*residualSet[fieldIndex]);

	soltransSet[fieldIndex]->prepare_for_coarsening_and_refinement(*solutionSet[fieldIndex]);
}
triangulation.execute_coarsening_and_refinement();

}

#include "../../include/matrixFreePDE_template_instantiations.h"
