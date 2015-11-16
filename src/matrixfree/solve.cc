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
  pcout << "\nsolving...\n\n";
  
  //time dependent BVP
  if (isTimeDependentBVP){
    //output initial conditions for time dependent BVP
    if (writeOutput) outputResults();
    
    //time stepping
    pcout << "\nTime stepping parameters: timeStep: " << dtValue << "  timeFinal: " << finalTime << "  timeIncrements: " << totalIncrements << "\n";

    for (currentIncrement=1; currentIncrement<=totalIncrements; ++currentIncrement){
      //increment current time
      currentTime+=dtValue;
      pcout << "\ntime increment:" << currentIncrement << "  time: " << currentTime << "\n";
      if (currentTime>=finalTime){
	pcout << "\ncurrentTime>=timeFinal. Ending time stepping\n";
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
    if (totalIncrements>1){
      pcout << "solve.h: this problem has only ELLIPTIC fields, hence neglecting totalIncrementsV>1 \n";
    }
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
