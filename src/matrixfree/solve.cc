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
    if (writeOutput) {
      outputResults();
#ifdef calc_energy
      if (calc_energy == true){
	computeEnergy();
	outputFreeEnergy(freeEnergyValues);
      }
#endif
    }
    
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
      
      //check and perform adaptive mesh refinement
      computing_timer.enter_section("matrixFreePDE: AMR");
      adaptiveRefine(currentIncrement);
      computing_timer.exit_section("matrixFreePDE: AMR");
 
      //solve time increment
      solveIncrement();

      //output results to file
      if ((writeOutput) && (currentIncrement%skipOutputSteps==0)){
    	  outputResults();
#ifdef calc_energy
	  if (calc_energy == true){
	    computeEnergy();
	    outputFreeEnergy(freeEnergyValues);
	  }
#endif
	  
      }
    }
  }
  //time independent BVP
  else{
    if (totalIncrements>1){
      pcout << "solve.h: this problem has only ELLIPTIC fields, hence neglecting totalIncrementsV>1 \n";
    }
    totalIncrements=1;

    //check and perform adaptive mesh refinement
    computing_timer.enter_section("matrixFreePDE: AMR");
    adaptiveRefine(0);
    computing_timer.exit_section("matrixFreePDE: AMR");
    
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
