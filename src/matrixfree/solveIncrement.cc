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

  //updateRHS
  updateRHS();
  
  //solve for each field
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    
    //Parabolic (first order derivatives in time) fields
    if (fields[fieldIndex].pdetype==PARABOLIC){
      //explicit-time step each DOF
      for (unsigned int dof=0; dof<solutionSet[fieldIndex]->local_size(); ++dof){
	solutionSet[fieldIndex]->local_element(dof)+=\
	  invM.local_element(dof)*timeStep*residualSet[fieldIndex]->local_element(dof);
      }
      char buffer[200];
      sprintf(buffer, "field '%s' [explicit solve]: current solution: %12.6e, current residual:%12.6e\n", \
	      fields[fieldIndex].name.c_str(),				\
	      solutionSet[fieldIndex]->l2_norm(),			\
	      residualSet[fieldIndex]->l2_norm()); 
      pcout<<buffer; 
    }
    
    //Elliptic (time-independent) fields
    else if (fields[fieldIndex].pdetype==ELLIPTIC){
      //implicit solve
#ifdef solverType
      SolverControl solver_control(maxSolverIterations, relSolverTolerance*residualSet[fieldIndex]->l2_norm());
      solverType<vectorType> solver(solver_control);
      //set implicitFieldIndex to mark field which is being implicitly
      //solved (required in the computeLHS() method)
      implicitFieldIndex=fieldIndex;
      try{
	solver.solve(*this, *solutionSet[fieldIndex], *residualSet[fieldIndex], IdentityMatrix(solutionSet[fieldIndex]->size()));
      }
      catch (...) {
	pcout << "\nWarning: solver did not converge as per set tolerances. consider increasing maxSolverIterations or decreasing relSolverTolerance.\n";
      }
      char buffer[200];
      sprintf(buffer, "field '%s' [implicit solve]: initial residual:%12.6e, current residual:%12.6e, nsteps:%u, tolerance criterion:%12.6e\n",\
	      fields[fieldIndex].name.c_str(),				\
	      solver_control.initial_value(),				\
	      solver_control.last_value(),				\
	      solver_control.last_step(), solver_control.tolerance()); 
      pcout<<buffer; 
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
  }
  pcout << "wall time: " << time.wall_time() << "s\n";
  //log time
  computing_timer.exit_section("matrixFreePDE: solveIncrements"); 
}

#endif
