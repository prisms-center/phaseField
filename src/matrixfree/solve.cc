//solve() method for MatrixFreePDE class

#ifndef SOLVE_MATRIXFREE_H
#define SOLVE_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//solve BVP
template <int dim>
void MatrixFreePDE<dim>::solve(){
  //log time
  computing_timer.enter_section("matrixFreePDE: solve"); 
  
  //time dependent BVP
  if (isTimeDependentBVP){
    //initialize time step variables
    timeStep=timeStepV;
    finalTime=finalTimeV;
    totalIncrements=totalIncrementsV;    
    //output initial conditions for time dependent BVP
    if (writeOutput) outputResults();

    //time step
    for (currentIncrement=1; currentIncrement<totalIncrements; ++currentIncrement){
      //increment current time
      currentTime+=timeStep;
      pcout << "\ntime increment:" << currentIncrement << "  time: " << currentTime << "\n";
      if (currentTime>=finalTime){
	pcout << "\ncurrentTime>=finalTime. Ending time stepping\n";
	break;
      }
      //solve time increment
      solveIncrement();
      //output results to file
      if ((writeOutput) && (currentIncrement%skipOutputSteps==0)){
	outputResults();
      }
    }
  }
  //time independent BVP
  else{
    totalIncrements=1;  
    //solve
    solveIncrement();
    //output results to file
    if ((writeOutput) && (currentIncrement%skipOutputSteps==0)){
      outputResults();
    }
  }

  //log time
  computing_timer.exit_section("matrixFreePDE: solve"); 
}

#endif
