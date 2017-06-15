//solveIncrement() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

//solve each time increment
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::solveIncrement(){
  //log time
  computing_timer.enter_section("matrixFreePDE: solveIncrements");
  Timer time;
  char buffer[200];

  // Get the list of new nuclei, if relevant
  // if (userInputs.nucleation_occurs){
  //  computing_timer.enter_section("matrixFreePDE: nucleation");
  //  getNucleiList();
  //  computing_timer.exit_section("matrixFreePDE: nucleation");
  // }

  //compute residual vectors
  computeRHS();

  //solve for each field
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
	  currentFieldIndex = fieldIndex; // Used in computeLHS()

    //Parabolic (first order derivatives in time) fields
    if (fields[fieldIndex].pdetype==PARABOLIC){

    	// Explicit-time step each DOF
    	// Takes advantage of knowledge that the length of solutionSet and residualSet is an integer multiple of the length of invM for vector variables
    	unsigned int invM_size = invM.local_size();

    	for (unsigned int dof=0; dof<solutionSet[fieldIndex]->local_size(); ++dof){

    		solutionSet[fieldIndex]->local_element(dof)=			\
    				invM.local_element(dof%invM_size)*residualSet[fieldIndex]->local_element(dof);
    	}

    	//apply constraints
    	constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
    	constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
    	//sync ghost DOF's
    	solutionSet[fieldIndex]->update_ghost_values();
    	//
    	if (currentIncrement%userInputs.skip_print_steps==0){
    		sprintf(buffer, "field '%2s' [explicit solve]: current solution: %12.6e, current residual:%12.6e\n", \
    				fields[fieldIndex].name.c_str(),				\
					solutionSet[fieldIndex]->l2_norm(),			\
					residualSet[fieldIndex]->l2_norm());
    		pcout<<buffer;
    	}
    }
    //Elliptic (time-independent) fields
    else if (fields[fieldIndex].pdetype==ELLIPTIC){

    	//implicit solve
		#ifdef solverType
		if (currentIncrement%skipImplicitSolves==0){
			//apply Dirichlet BC's
			// Loops through all DoF to which ones have Dirichlet BCs applied, replace the ones that do with the Dirichlet value
			// Is this needed? Why are we applying BCs to the residualSet?
			// This clears the residual where we want to apply Dirichlet BCs, otherwise the solver sees a positive residual
			for (std::map<types::global_dof_index, double>::const_iterator it=valuesDirichletSet[fieldIndex]->begin(); it!=valuesDirichletSet[fieldIndex]->end(); ++it){
				if (residualSet[fieldIndex]->in_local_range(it->first)){
					(*residualSet[fieldIndex])(it->first) = 0.0; //it->second; //*jacobianDiagonal(it->first);
				}
			}

			//solver controls
			double tol_value;
			if (userInputs.abs_tol == true){
				tol_value = userInputs.solver_tolerance;
				//if (currentIncrement > 1) tol_value = 1000000.0;
			}
			else {
				tol_value = userInputs.solver_tolerance*residualSet[fieldIndex]->l2_norm();
			}

			SolverControl solver_control(userInputs.max_solver_iterations, tol_value);

			// Currently the only allowed solver is SolverCG, the SolverType input variable is a dummy
			SolverCG<vectorType> solver(solver_control);

			//solve
			try{
				if (fields[fieldIndex].type == SCALAR){
					dU_scalar=0.0;
					solver.solve(*this, dU_scalar, *residualSet[fieldIndex], IdentityMatrix(solutionSet[fieldIndex]->size()));
				}
				else {
					dU_vector=0.0;
					solver.solve(*this, dU_vector, *residualSet[fieldIndex], IdentityMatrix(solutionSet[fieldIndex]->size()));
				}
			}
			catch (...) {
				pcout << "\nWarning: implicit solver did not converge as per set tolerances. consider increasing maxSolverIterations or decreasing solverTolerance.\n";
			}
			if (fields[fieldIndex].type == SCALAR){
				*solutionSet[fieldIndex]+=dU_scalar;
			}
			else {
				*solutionSet[fieldIndex]+=dU_vector;
			}

			// Apply hanging node and periodic constraints
			constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
			//sync ghost DOF's
			solutionSet[fieldIndex]->update_ghost_values();
			//
			 if (currentIncrement%userInputs.skip_print_steps==0){
				 double dU_norm;
				 if (fields[fieldIndex].type == SCALAR){
					 dU_norm = dU_scalar.l2_norm();
				 }
				 else {
					 dU_norm = dU_vector.l2_norm();
				 }
				 sprintf(buffer, "field '%2s' [implicit solve]: initial residual:%12.6e, current residual:%12.6e, nsteps:%u, tolerance criterion:%12.6e, solution: %12.6e, dU: %12.6e\n", \
						 fields[fieldIndex].name.c_str(),			\
						 residualSet[fieldIndex]->l2_norm(),			\
						 solver_control.last_value(),				\
						 solver_control.last_step(), solver_control.tolerance(), solutionSet[fieldIndex]->l2_norm(), dU_norm);
				 pcout<<buffer;
			 }
		}
		else{
			sprintf(buffer, "field '%2s' [implicit solve]: current residual:%12.6e\n", \
					fields[fieldIndex].name.c_str(),			\
					residualSet[fieldIndex]->l2_norm());
			pcout<<buffer;
			pcout << "skipping implicit solve as currentIncrement%skipImplicitSolves!=0\n";
		}

		#else
		pcout << "\nError: solverType not defined. This is required for ELLIPTIC fields.\n\n";
		exit (-1);
		#endif
	  }

	  //Hyperbolic (second order derivatives in time) fields and general
	  //non-linear PDE types not yet implemented
	  else{
		  pcout << "matrixFreePDE.h: unknown field pdetype\n";
		  exit(-1);
	  }

	  //check if solution is nan
	  if (!numbers::is_finite(solutionSet[fieldIndex]->l2_norm())){
		  sprintf(buffer, "ERROR: field '%s' solution is NAN. exiting.\n\n",
				  fields[fieldIndex].name.c_str());
		  pcout<<buffer;
		  exit(-1);
	  }
  }
  if (currentIncrement%userInputs.skip_print_steps==0){
  pcout << "wall time: " << time.wall_time() << "s\n";
  }
  //log time
  computing_timer.exit_section("matrixFreePDE: solveIncrements");
}

#include "../../include/matrixFreePDE_template_instantiations.h"
