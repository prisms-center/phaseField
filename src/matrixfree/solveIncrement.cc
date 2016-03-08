//solveIncrement() method for MatrixFreePDE class

#ifndef SOLVEINCREMENT_MATRIXFREE_H
#define SOLVEINCREMENT_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//solve each time increment
template <int dim>
void MatrixFreePDE<dim>::solveIncrement(){
  //log time
  computing_timer.enter_section("matrixFreePDE: solveIncrements");
  Timer time; 
  char buffer[200];

  //modify fields (rarely used. Typically used in problems involving nucleation)
  if (nucleation_occurs == true) modifySolutionFields();

  //compute residual vectors
  computeRHS();

  //solve for each field
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
	  //Parabolic (first order derivatives in time) fields
	  if (fields[fieldIndex].pdetype==PARABOLIC){
		  //explicit-time step each DOF
		  for (unsigned int dof=0; dof<solutionSet[fieldIndex]->local_size(); ++dof){
			  solutionSet[fieldIndex]->local_element(dof)=			\
					  invM.local_element(dof)*residualSet[fieldIndex]->local_element(dof);
		  }
//      sprintf(buffer, "field '%2s' [explicit solve]: current residual:%12.6e\n", \
//	      fields[fieldIndex].name.c_str(),				\
//	      residualSet[fieldIndex]->l2_norm());
//      pcout<<buffer;
	  }
	  //Elliptic (time-independent) fields
	  else if (fields[fieldIndex].pdetype==ELLIPTIC){
		  //implicit solve
#ifdef solverType
		  #if abs_tol == true
			  SolverControl solver_control(maxSolverIterations, absSolverTolerance);

		  #else
			  SolverControl solver_control(maxSolverIterations, relSolverTolerance*residualSet[fieldIndex]->l2_norm());
		  #endif
		  solverType<vectorType> solver(solver_control);
		  if (currentIncrement%skipImplicitSolves==0){
			  try{
				  solver.solve(*this, dU, *residualSet[fieldIndex], IdentityMatrix(solutionSet[fieldIndex]->size()));
			  }
			  catch (...) {
				  pcout << "\nWarning: implicit solver did not converge as per set tolerances. consider increasing maxSolverIterations or decreasing relSolverTolerance.\n";
			  }
			  *solutionSet[fieldIndex]+=dU;
			  sprintf(buffer, "field '%2s' [implicit solve]: initial residual:%12.6e, current residual:%12.6e, nsteps:%u, tolerance criterion:%12.6e\n", \
					  fields[fieldIndex].name.c_str(),			\
					  residualSet[fieldIndex]->l2_norm(),			\
					  solver_control.last_value(),				\
					  solver_control.last_step(), solver_control.tolerance());
			  pcout<<buffer;
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
	  if (!numbers::is_finite(solutionSet[fieldIndex]->norm_sqr())){
		  sprintf(buffer, "ERROR: field '%s' solution is NAN. exiting.\n\n",
				  fields[fieldIndex].name.c_str());
		  pcout<<buffer;
		  exit(-1);
	  }
  }
  pcout << "wall time: " << time.wall_time() << "s\n";
  //log time
  computing_timer.exit_section("matrixFreePDE: solveIncrements"); 
}

#endif
