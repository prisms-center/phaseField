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

  getOutputTimeSteps(outputCondition,numOutputs,userInputs.user_given_time_step_list,outputTimeStepList);
  int currentOutput = 0;

  //time dependent BVP
  if (isTimeDependentBVP){
    //output initial conditions for time dependent BVP
	  if ((writeOutput) && (outputTimeStepList[currentOutput] == 0)) {

			  outputResults();
			  if (userInputs.calc_energy == true){
				  computeEnergy();
				  outputFreeEnergy(freeEnergyValues);
			  }
			  currentOutput++;
    }
    
    //time stepping
    pcout << "\nTime stepping parameters: timeStep: " << userInputs.dtValue << "  timeFinal: " << userInputs.finalTime << "  timeIncrements: " << userInputs.totalIncrements << "\n";
    
    for (currentIncrement=1; currentIncrement<=userInputs.totalIncrements; ++currentIncrement){
      //increment current time
      currentTime+=userInputs.dtValue;
      if (currentIncrement%skipPrintSteps==0){
      pcout << "\ntime increment:" << currentIncrement << "  time: " << currentTime << "\n";
      }

      //check and perform adaptive mesh refinement
      computing_timer.enter_section("matrixFreePDE: AMR");
      adaptiveRefine(currentIncrement);
      computing_timer.exit_section("matrixFreePDE: AMR");
 
      //solve time increment
      solveIncrement();

      //Ghost the solution vectors. Also apply the Dirichet BC's (if any) on the solution vectors
      for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    	  constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
    	  solutionSet[fieldIndex]->update_ghost_values();
      }

      //output results to file
      if ((writeOutput) && (outputTimeStepList[currentOutput] == currentIncrement)) {
    	  outputResults();

    	  if (userInputs.calc_energy == true){
    		  computeEnergy();
    		  outputFreeEnergy(freeEnergyValues);
    	  }

    	  currentOutput++;
      }
    }
  }

  //time independent BVP
  else{
    if (userInputs.totalIncrements>1){
      pcout << "solve.h: this problem has only ELLIPTIC fields, hence neglecting totalIncrementsV>1 \n";
    }
    userInputs.totalIncrements=1;

    //check and perform adaptive mesh refinement
    computing_timer.enter_section("matrixFreePDE: AMR");
    //adaptiveRefine(0);
    computing_timer.exit_section("matrixFreePDE: AMR");
    
    //solve
    solveIncrement();

    //Ghost the solution vectors. Also apply the Dirichet BC's (if any) on the solution vectors
    for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    	constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
    	solutionSet[fieldIndex]->update_ghost_values();
    }

    //output results to file
    if (writeOutput){
    	outputResults();
    	if (userInputs.calc_energy == true){
    		computeEnergy();
    		outputFreeEnergy(freeEnergyValues);
    	}
    }

  }

  //log time
  computing_timer.exit_section("matrixFreePDE: solve"); 
}

#endif
