//solveIncrement() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

//solve each time increment
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::solveIncrement(){

    //log time
    computing_timer.enter_section("matrixFreePDE: solveIncrements");
    Timer time;
    char buffer[200];

    //compute residual vectors
    computeRHS();

    //solve for each field
    for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
        currentFieldIndex = fieldIndex; // Used in computeLHS()

        // Add Neumann BC terms to the residual vector for the current field, if appropriate
        // Currently commented out because it isn't working yet
        //applyNeumannBCs();

        //Parabolic (first order derivatives in time) fields
        if (fields[fieldIndex].pdetype==PARABOLIC){

            // Explicit-time step each DOF
            // Takes advantage of knowledge that the length of solutionSet and residualSet is an integer multiple of the length of invM for vector variables
            unsigned int invM_size = invM.local_size();
            for (unsigned int dof=0; dof<solutionSet[fieldIndex]->local_size(); ++dof){
                solutionSet[fieldIndex]->local_element(dof)=			\
                invM.local_element(dof%invM_size)*residualSet[fieldIndex]->local_element(dof);
            }

            // Set the Dirichelet values (hanging node constraints don't need to be distributed every time step, only at output)
            constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
            solutionSet[fieldIndex]->update_ghost_values();

            // Print update to screen
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
            unsigned int allowed_solver_iterations;
            double residual_norm = residualSet[fieldIndex]->l2_norm();
            if (userInputs.abs_tol == true){

                if (residual_norm > userInputs.solver_tolerance){
                    tol_value = userInputs.solver_tolerance;
                    allowed_solver_iterations = userInputs.max_solver_iterations;
                }
                else {
                    tol_value = 0.0;
                    allowed_solver_iterations = 1;
                }

                tol_value = userInputs.solver_tolerance;
            }
            else {
                tol_value = userInputs.solver_tolerance*residual_norm;
            }

            IterationNumberControl solver_control(allowed_solver_iterations, tol_value);

            // Currently the only allowed solver is SolverCG, the SolverType input variable is a dummy
            SolverCG<vectorType> solver(solver_control);

            //solve

            if (fields[fieldIndex].type == SCALAR){
                dU_scalar=0.0;
                solver.solve(*this, dU_scalar, *residualSet[fieldIndex], IdentityMatrix(solutionSet[fieldIndex]->size()));
            }
            else {
                dU_vector=0.0;
                solver.solve(*this, dU_vector, *residualSet[fieldIndex], IdentityMatrix(solutionSet[fieldIndex]->size()));
            }

            if (solver_control.last_step() == userInputs.max_solver_iterations){
                pcout << "\nWarning: Implicit solver did not converge as per set tolerances. Consider increasing maxSolverIterations or decreasing solverTolerance. " << solver_control.last_step()<< "steps were taken.\n";
            }

            if (fields[fieldIndex].type == SCALAR){
                *solutionSet[fieldIndex]+=dU_scalar;
            }
            else {
                *solutionSet[fieldIndex]+=dU_vector;
            }

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
                residual_norm,			\
                solver_control.last_value(),				\
                solver_control.last_step(), solver_control.tolerance(), solutionSet[fieldIndex]->l2_norm(), dU_norm);
                pcout<<buffer;
            }

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
