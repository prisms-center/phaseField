//mesh refinement methods for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

//default implementation of adaptive mesh refinement 
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::adaptiveRefine(unsigned int currentIncrement){
if (userInputs.h_adaptivity == true){
	if ( (currentIncrement == 0) ){
		for (unsigned int remesh_index=0; remesh_index < (userInputs.max_refinement_level-userInputs.min_refinement_level); remesh_index++){
			reinit();
		}
	}
	else if ( (currentIncrement%userInputs.skip_remeshing_steps==0) ){
		reinit();
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


QGauss<dim>  quadrature(degree+1);
FEValues<dim> fe_values (*FESet[userInputs.refine_criterion_fields[0]], quadrature, update_values);
const unsigned int   num_quad_points = quadrature.size();

std::vector<double> errorOut(num_quad_points);

typename DoFHandler<dim>::active_cell_iterator cell = dofHandlersSet_nonconst[userInputs.refine_criterion_fields[0]]->begin_active(), endc = dofHandlersSet_nonconst[userInputs.refine_criterion_fields[0]]->end();

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

		if ( (mark_refine == true) ){
			cell->set_refine_flag();
		}
		else {
			cell->set_coarsen_flag();
		}
	}
}

}


//refine grid method
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::refineGrid (){

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
  if ((unsigned int)cell->level() <= userInputs.min_refinement_level ){
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

}

#include "../../include/matrixFreePDE_template_instantiations.h"

