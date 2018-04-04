//solve() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

//solve BVP
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::solve(){
    //log time
    computing_timer.enter_section("matrixFreePDE: solve");
    pcout << "\nsolving...\n\n";

    //time dependent BVP
    if (isTimeDependentBVP){

        //time stepping
        pcout << "\nTime stepping parameters: timeStep: " << userInputs.dtValue << "  timeFinal: " << userInputs.finalTime << "  timeIncrements: " << userInputs.totalIncrements << "\n";

        // Do an initial solve to set the elliptic fields
        solveIncrement(true);

        //output initial conditions for time dependent BVP
        if (userInputs.outputTimeStepList[currentOutput] == currentIncrement) {

            for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
                constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                solutionSet[fieldIndex]->update_ghost_values();
            }
            outputResults();
            currentOutput++;
        }

        if (userInputs.checkpointTimeStepList[currentCheckpoint] == currentIncrement) {
            save_checkpoint();
            currentCheckpoint++;
        }

        // Increase the current increment from 0 to 1 now that the initial conditions have been output
        currentIncrement++;

        // Cycle up to the proper output and checkpoint counters
        while (userInputs.outputTimeStepList.size() > 0 && userInputs.outputTimeStepList[currentOutput] < currentIncrement){
            currentOutput++;
        }
        while (userInputs.checkpointTimeStepList.size() > 0 && userInputs.checkpointTimeStepList[currentCheckpoint] < currentIncrement){
            currentCheckpoint++;
        }

        // This is the main time-stepping loop
        for (; currentIncrement<=userInputs.totalIncrements; ++currentIncrement){
            //increment current time
            currentTime+=userInputs.dtValue;
            if (currentIncrement%userInputs.skip_print_steps==0){
                pcout << "\ntime increment:" << currentIncrement << "  time: " << currentTime << "\n";
            }

            //check and perform adaptive mesh refinement
            adaptiveRefine(currentIncrement);

            // Update the list of nuclei (if relevant)
            updateNucleiList();

            //solve time increment
            solveIncrement(false);

            // Output results to file (on the proper increments)
            if (userInputs.outputTimeStepList[currentOutput] == currentIncrement) {
                for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
                    constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                    constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                    solutionSet[fieldIndex]->update_ghost_values();
                }
                outputResults();
                currentOutput++;
            }

            // Create a checkpoint (on the proper increments)
            if (userInputs.checkpointTimeStepList[currentCheckpoint] == currentIncrement) {
                save_checkpoint();
                currentCheckpoint++;
            }

        }
    }

    //time independent BVP
    else{
        //solve
        solveIncrement(false);

        //output results to file
        outputResults();
    }

    //log time
    computing_timer.exit_section("matrixFreePDE: solve");
}

#include "../../include/matrixFreePDE_template_instantiations.h"
