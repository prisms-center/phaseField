// solve() method for MatrixFreePDE class

#include <core/matrixFreePDE.h>

// solve BVP
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::solve()
{
  // log time
  computing_timer.enter_subsection("matrixFreePDE: solve");
  pcout << "\nsolving...\n\n";

  bool nucleation_complete = false;

  // time dependent BVP
  if (isTimeDependentBVP)
    {
      // If grain reassignment is activated, reassign grains
      if (userInputs.grain_remapping_activated and
          (currentIncrement % userInputs.skip_grain_reassignment_steps == 0 or
           currentIncrement == 0))
        {
          reassignGrains();
        }

      // For any nonlinear equation, set the initial guess as the solution to
      // Laplace's equations
      generatingInitialGuess = true;
      setNonlinearEqInitialGuess();
      generatingInitialGuess = false;

      // Do an initial solve to set the elliptic fields
      solveIncrement(true);

      // output initial conditions for time dependent BVP
      if (userInputs.outputTimeStepList[currentOutput] == currentIncrement)
        {
          for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
            {
              constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
              constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
              solutionSet[fieldIndex]->update_ghost_values();
            }
          outputResults();
          currentOutput++;
        }

      if (userInputs.checkpointTimeStepList[currentCheckpoint] == currentIncrement)
        {
          save_checkpoint();
          currentCheckpoint++;
        }

      // Increase the current increment from 0 to 1 now that the initial
      // conditions have been output
      currentIncrement++;

      // Cycle up to the proper output and checkpoint counters
      while (userInputs.outputTimeStepList.size() > 0 &&
             userInputs.outputTimeStepList[currentOutput] < currentIncrement)
        {
          currentOutput++;
        }
      while (userInputs.checkpointTimeStepList.size() > 0 &&
             userInputs.checkpointTimeStepList[currentCheckpoint] < currentIncrement)
        {
          currentCheckpoint++;
        }

      // time stepping
      pcout << "\nTime stepping parameters: timeStep: " << userInputs.dtValue
            << "  timeFinal: " << userInputs.finalTime
            << "  timeIncrements: " << userInputs.totalIncrements << "\n";

      // This is the main time-stepping loop
      for (; currentIncrement <= userInputs.totalIncrements; ++currentIncrement)
        {
          // increment current time
          currentTime += userInputs.dtValue;
          if (currentIncrement % userInputs.skip_print_steps == 0)
            {
              pcout << "\ntime increment:" << currentIncrement
                    << "  time: " << currentTime << "\n";
            }

          // check and perform adaptive mesh refinement
          if (userInputs.h_adaptivity == true &&
              currentIncrement % userInputs.skip_remeshing_steps == 0)
            {
              computing_timer.enter_subsection("matrixFreePDE: AMR");

              AMR.do_adaptive_refinement(currentIncrement);
              reinit();

              computing_timer.leave_subsection("matrixFreePDE: AMR");
            }

          // Update the list of nuclei (if relevant)
          updateNucleiList();

          if (currentTime >= userInputs.nucleation_end_time && !nucleation_complete)
            {
              pcout << "Total number of nuclei = " << nuclei.size() << std::endl;
	      nucleation_complete = true;
            }

          // If grain reassignment is activated, reassign grains
          if (userInputs.grain_remapping_activated and
              (currentIncrement % userInputs.skip_grain_reassignment_steps == 0 or
               currentIncrement == 0))
            {
              reassignGrains();
            }

          // solve time increment
          solveIncrement(false);

          // Output results to file (on the proper increments)
          if (userInputs.outputTimeStepList[currentOutput] == currentIncrement)
            {
              for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
                {
                  constraintsDirichletSet[fieldIndex]->distribute(
                    *solutionSet[fieldIndex]);
                  constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                  solutionSet[fieldIndex]->update_ghost_values();
                }
              outputResults();
              if (userInputs.print_timing_with_output &&
                  currentIncrement < userInputs.totalIncrements)
                {
                  computing_timer.print_summary();
                }

              currentOutput++;
            }

          // Create a checkpoint (on the proper increments)
          if (userInputs.checkpointTimeStepList[currentCheckpoint] == currentIncrement)
            {
              save_checkpoint();
              currentCheckpoint++;
            }
        }
    }

  // time independent BVP
  else
    {
      generatingInitialGuess = false;

      // solve
      solveIncrement(false);

      // output results to file
      outputResults();
    }

  // log time
  computing_timer.leave_subsection("matrixFreePDE: solve");
}

template class MatrixFreePDE<2, 1>;
template class MatrixFreePDE<3, 1>;

template class MatrixFreePDE<2, 2>;
template class MatrixFreePDE<3, 2>;

template class MatrixFreePDE<3, 3>;
template class MatrixFreePDE<2, 3>;

template class MatrixFreePDE<3, 4>;
template class MatrixFreePDE<2, 4>;

template class MatrixFreePDE<3, 5>;
template class MatrixFreePDE<2, 5>;

template class MatrixFreePDE<3, 6>;
template class MatrixFreePDE<2, 6>;
