//solve() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

//solve BVP
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::solve(){
  //log time
  computing_timer.enter_section("matrixFreePDE: solve");
  pcout << "\nsolving...\n\n";

  getOutputTimeSteps(userInputs.output_condition,userInputs.num_outputs,userInputs.user_given_time_step_list,outputTimeStepList);
  int currentOutput = 0;

  //time dependent BVP
  if (isTimeDependentBVP){
    //output initial conditions for time dependent BVP
	  if ((writeOutput) && (outputTimeStepList[currentOutput] == 0)) {

          for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
              constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
              constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
              solutionSet[fieldIndex]->update_ghost_values();
          }
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
      if (currentIncrement%userInputs.skip_print_steps==0){
      pcout << "\ntime increment:" << currentIncrement << "  time: " << currentTime << "\n";
      }

      //check and perform adaptive mesh refinement
      computing_timer.enter_section("matrixFreePDE: AMR");
      adaptiveRefine(currentIncrement);
      computing_timer.exit_section("matrixFreePDE: AMR");

      //solve time increment
      solveIncrement();

      //output results to file
      if ((writeOutput) && (outputTimeStepList[currentOutput] == currentIncrement)) {
          for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
              constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
              constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
              solutionSet[fieldIndex]->update_ghost_values();
          }
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
    //computing_timer.enter_section("matrixFreePDE: AMR");
    //adaptiveRefine(0);
    //computing_timer.exit_section("matrixFreePDE: AMR");

    //solve
    solveIncrement();

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

#include "../../include/matrixFreePDE_template_instantiations.h"
