//solveIncrement() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

//solve each time increment
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::solveIncrement(bool skip_time_dependent){

    //log time
    computing_timer.enter_section("matrixFreePDE: solveIncrements");
    Timer time;
    char buffer[200];

    // Get the RHS of the equations
    // Ideally this would be just for the explicit fields, but for now this is all of them
    computeRHS();

    //solve for each field
    for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
        currentFieldIndex = fieldIndex; // Used in computeLHS()

        // Add Neumann BC terms to the residual vector for the current field, if appropriate
        // Currently commented out because it isn't working yet
        //applyNeumannBCs();

        //Parabolic (first order derivatives in time) fields
        if (fields[fieldIndex].pdetype==PARABOLIC && !skip_time_dependent){

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

        //check if solution is nan
        if (!numbers::is_finite(solutionSet[fieldIndex]->l2_norm())){
            sprintf(buffer, "ERROR: field '%s' solution is NAN. exiting.\n\n",
            fields[fieldIndex].name.c_str());
            pcout<<buffer;
            exit(-1);
        }
    }

    // Now, update the non-explicit variables
    // For the time being, this is just the elliptic equations, but implicit parabolic and auxilary equations should also be here
    if (isEllipticBVP){
        double nonlinear_convergence_tolerance = 1.0e-4; // For testing purposes, this is hardcoded
        unsigned int max_nonlinear_it = 10;
        bool nonlinear_it_converged = false;
        unsigned int nonlinear_it_index = 0;

        while (!nonlinear_it_converged){
            nonlinear_it_converged = true; // Set to true here and will be set to false if any variable isn't converged

            // Update residualSet for the non-explicitly updated variables
            //compute_nonexplicit_RHS()
            // Ideally, I'd just do this for the non-explicit variables, but for now I'll do all of them
            // this is a little redundant, but hopefully not too terrible
            computeRHS();

            for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
                currentFieldIndex = fieldIndex; // Used in computeLHS()

                if (fields[fieldIndex].pdetype==ELLIPTIC){
                    dealii::parallel::distributed::Vector<double> solution_diff = *solutionSet[fieldIndex];

                    //apply Dirichlet BC's
                    // Loops through all DoF to which ones have Dirichlet BCs applied, replace the ones that do with the Dirichlet value
                    // This clears the residual where we want to apply Dirichlet BCs, otherwise the solver sees a positive residual
                    for (std::map<types::global_dof_index, double>::const_iterator it=valuesDirichletSet[fieldIndex]->begin(); it!=valuesDirichletSet[fieldIndex]->end(); ++it){
                        if (residualSet[fieldIndex]->in_local_range(it->first)){
                            (*residualSet[fieldIndex])(it->first) = 0.0;
                        }
                    }

                    //solver controls
                    double tol_value;
                    if (userInputs.abs_tol == true){
                        tol_value = userInputs.solver_tolerance;
                    }
                    else {
                        tol_value = userInputs.solver_tolerance*residualSet[fieldIndex]->l2_norm();
                    }

                    SolverControl solver_control(userInputs.max_solver_iterations, tol_value);

                    // Currently the only allowed solver is SolverCG, the SolverType input variable is a dummy
                    SolverCG<vectorType> solver(solver_control);

                    //solve
                    try{
                        if (fields[fieldIndex].type == SCALAR){
                            dU_scalar=0.0;
                            solver.solve(*this, dU_scalar, *residualSet[fieldIndex], IdentityMatrix(solutionSet[fieldIndex]->size()));
                        }
                        else {
                            dU_vector=0.0;
                            solver.solve(*this, dU_vector, *residualSet[fieldIndex], IdentityMatrix(solutionSet[fieldIndex]->size()));
                        }
                    }
                    catch (...) {
                        pcout << "\nWarning: implicit solver did not converge as per set tolerances. consider increasing maxSolverIterations or decreasing solverTolerance.\n";
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
                        residualSet[fieldIndex]->l2_norm(),			\
                        solver_control.last_value(),				\
                        solver_control.last_step(), solver_control.tolerance(), solutionSet[fieldIndex]->l2_norm(), dU_norm);
                        pcout<<buffer;
                    }

                    // Check to see if this individual variable has converged
                    double diff;
                    if (fields[fieldIndex].type == SCALAR){
                        diff = dU_scalar.linfty_norm();
                    }
                    else {
                        diff = dU_vector.linfty_norm();
                    }
                    pcout << "Max difference between nonlinear iterations: " << diff << " " << nonlinear_it_index << " " << currentIncrement << std::endl;

                    if (diff > nonlinear_convergence_tolerance && nonlinear_it_index < max_nonlinear_it){
                        nonlinear_it_converged = false;
                    }
                }
            }

            nonlinear_it_index++;
        }
    }

    if (currentIncrement%userInputs.skip_print_steps==0){
        pcout << "wall time: " << time.wall_time() << "s\n";
    }
    //log time
    computing_timer.exit_section("matrixFreePDE: solveIncrements");

}

#include "../../include/matrixFreePDE_template_instantiations.h"
