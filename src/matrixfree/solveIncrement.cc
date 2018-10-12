//solveIncrement() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"
#include <deal.II/lac/solver_cg.h>

//solve each time increment
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::solveIncrement(bool skip_time_dependent){

    //log time
    computing_timer.enter_section("matrixFreePDE: solveIncrements");
    Timer time;
    char buffer[200];

    // Get the RHS of the explicit equations
    if (hasExplicitEquation && !skip_time_dependent){
        computeExplicitRHS();
    }


    //solve for each field
    for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
        currentFieldIndex = fieldIndex; // Used in computeLHS()

        // Add Neumann BC terms to the residual vector for the current field, if appropriate
        // Currently commented out because it isn't working yet
        //applyNeumannBCs();

        //Parabolic (first order derivatives in time) fields
        if (fields[fieldIndex].pdetype==EXPLICIT_TIME_DEPENDENT && !skip_time_dependent){

            // Explicit-time step each DOF
            // Takes advantage of knowledge that the length of solutionSet and residualSet is an integer multiple of the length of invM for vector variables
            unsigned int invM_size = invM.local_size();
            for (unsigned int dof=0; dof<solutionSet[fieldIndex]->local_size(); ++dof){
                solutionSet[fieldIndex]->local_element(dof)=			\
                invM.local_element(dof%invM_size)*residualSet[fieldIndex]->local_element(dof);
            }
            // Set the Dirichelet values (hanging node constraints don't need to be distributed every time step, only at output)
            if (has_Dirichlet_BCs){
                constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
            }
            //computing_timer.enter_section("matrixFreePDE: updateExplicitGhosts");
            solutionSet[fieldIndex]->update_ghost_values();
            //computing_timer.exit_section("matrixFreePDE: updateExplicitGhosts");

            // Print update to screen and confirm that solution isn't nan
            if (currentIncrement%userInputs.skip_print_steps==0){
                double solution_L2_norm = solutionSet[fieldIndex]->l2_norm();

                sprintf(buffer, "field '%2s' [explicit solve]: current solution: %12.6e, current residual:%12.6e\n", \
                fields[fieldIndex].name.c_str(),				\
                solution_L2_norm,			\
                residualSet[fieldIndex]->l2_norm());
                pcout<<buffer;

                if (!numbers::is_finite(solution_L2_norm)){
                    sprintf(buffer, "ERROR: field '%s' solution is NAN. exiting.\n\n",
                    fields[fieldIndex].name.c_str());
                    pcout<<buffer;
                    exit(-1);
               }

            }
        }

    }

    // Now, update the non-explicit variables
    // For the time being, this is just the elliptic equations, but implicit parabolic and auxilary equations should also be here
    if (hasNonExplicitEquation){

        bool nonlinear_it_converged = false;
        unsigned int nonlinear_it_index = 0;

        while (!nonlinear_it_converged){
            nonlinear_it_converged = true; // Set to true here and will be set to false if any variable isn't converged

            // Update residualSet for the non-explicitly updated variables
            //compute_nonexplicit_RHS()
            // Ideally, I'd just do this for the non-explicit variables, but for now I'll do all of them
            // this is a little redundant, but hopefully not too terrible
            computeNonexplicitRHS();

            for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
                currentFieldIndex = fieldIndex; // Used in computeLHS()

                if ( (fields[fieldIndex].pdetype == IMPLICIT_TIME_DEPENDENT && !skip_time_dependent) || fields[fieldIndex].pdetype == TIME_INDEPENDENT){

                    if (currentIncrement%userInputs.skip_print_steps==0 && userInputs.var_nonlinear[fieldIndex]){
                        sprintf(buffer, "field '%2s' [nonlinear solve]: current solution: %12.6e, current residual:%12.6e\n", \
                        fields[fieldIndex].name.c_str(),				\
                        solutionSet[fieldIndex]->l2_norm(),			\
                        residualSet[fieldIndex]->l2_norm());
                        pcout<<buffer;
                    }

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
                    if (userInputs.linear_solver_parameters.getToleranceType(fieldIndex) == ABSOLUTE_RESIDUAL){
                        tol_value = userInputs.linear_solver_parameters.getToleranceValue(fieldIndex);
                    }
                    else {
                        tol_value = userInputs.linear_solver_parameters.getToleranceValue(fieldIndex)*residualSet[fieldIndex]->l2_norm();
                    }

                    SolverControl solver_control(userInputs.linear_solver_parameters.getMaxIterations(fieldIndex), tol_value);

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

                    if (userInputs.var_nonlinear[fieldIndex]){

                        // Now that we have the calculated change in the solution, we need to select a damping coefficient
                        double damping_coefficient;

                        if (userInputs.nonlinear_solver_parameters.getBacktrackDampingFlag(fieldIndex)){
                            vectorType solutionSet_old = *solutionSet[fieldIndex];
                            double residual_old = residualSet[fieldIndex]->l2_norm();

                            damping_coefficient = 1.0;
                            bool damping_coefficient_found = false;
                            while (!damping_coefficient_found){
                                if (fields[fieldIndex].type == SCALAR){
                                    solutionSet[fieldIndex]->sadd(1.0,damping_coefficient,dU_scalar);
                                }
                                else {
                                    solutionSet[fieldIndex]->sadd(1.0,damping_coefficient,dU_vector);
                                }

                                computeNonexplicitRHS();

                                for (std::map<types::global_dof_index, double>::const_iterator it=valuesDirichletSet[fieldIndex]->begin(); it!=valuesDirichletSet[fieldIndex]->end(); ++it){
                                    if (residualSet[fieldIndex]->in_local_range(it->first)){
                                        (*residualSet[fieldIndex])(it->first) = 0.0;
                                    }
                                }

                                double residual_new = residualSet[fieldIndex]->l2_norm();

                                if (currentIncrement%userInputs.skip_print_steps==0){
                                    pcout << "    Old residual: " << residual_old << " Damping Coeff: " << damping_coefficient << " New Residual: " << residual_new << std::endl;
                                }

                                // An improved approach would use the Armijo–Goldstein condition to ensure a sufficent decrease in the residual. This way is just scales the residual.
                                if ( (residual_new < (residual_old*userInputs.nonlinear_solver_parameters.getBacktrackResidualDecreaseCoeff(fieldIndex))) || damping_coefficient < 1.0e-4){
                                    damping_coefficient_found = true;
                                }
                                else{
                                    damping_coefficient *= userInputs.nonlinear_solver_parameters.getBacktrackStepModifier(fieldIndex);
                                    *solutionSet[fieldIndex] = solutionSet_old;
                                }
                            }
                        }
                        else{
                            damping_coefficient = userInputs.nonlinear_solver_parameters.getDefaultDampingCoefficient(fieldIndex);

                            if (fields[fieldIndex].type == SCALAR){
                                solutionSet[fieldIndex]->sadd(1.0,damping_coefficient,dU_scalar);
                            }
                            else {
                                solutionSet[fieldIndex]->sadd(1.0,damping_coefficient,dU_vector);
                            }
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
                        if (userInputs.nonlinear_solver_parameters.getToleranceType(fieldIndex) == ABSOLUTE_SOLUTION_CHANGE){
                            double diff;

                            if (fields[fieldIndex].type == SCALAR){
                                diff = dU_scalar.l2_norm();
                            }
                            else {
                                diff = dU_vector.l2_norm();
                            }
                            if (currentIncrement%userInputs.skip_print_steps==0){
                                pcout << "Relative difference between nonlinear iterations: " << diff << " " << nonlinear_it_index << " " << currentIncrement << std::endl;
                            }

                            if (diff > userInputs.nonlinear_solver_parameters.getToleranceValue(fieldIndex) && nonlinear_it_index < userInputs.nonlinear_solver_parameters.getMaxIterations()){
                                nonlinear_it_converged = false;
                            }
                        }
                        else {
                            std::cerr << "PRISMS-PF Error: Nonlinear solver tolerance types other than ABSOLUTE_CHANGE have yet to be implemented." << std::endl;
                        }
                    }
                    else {
                        if (nonlinear_it_index ==0){

                            if (fields[fieldIndex].type == SCALAR){
                                *solutionSet[fieldIndex] += dU_scalar;
                            }
                            else {
                                *solutionSet[fieldIndex] += dU_vector;
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
                        }

                    }
                }
                else if (fields[fieldIndex].pdetype == AUXILIARY){

                    if (userInputs.var_nonlinear[fieldIndex] || nonlinear_it_index == 0){
                        // If the equation for this field is nonlinear, save the old solution
                        if (userInputs.var_nonlinear[fieldIndex]){
                            if (fields[fieldIndex].type == SCALAR){
                                dU_scalar = *solutionSet[fieldIndex];
                            }
                            else {
                                dU_vector = *solutionSet[fieldIndex];
                            }
                        }

                        // Explicit-time step each DOF
                        // Takes advantage of knowledge that the length of solutionSet and residualSet is an integer multiple of the length of invM for vector variables
                        unsigned int invM_size = invM.local_size();
                        for (unsigned int dof=0; dof<solutionSet[fieldIndex]->local_size(); ++dof){
                            solutionSet[fieldIndex]->local_element(dof)=			\
                            invM.local_element(dof%invM_size)*residualSet[fieldIndex]->local_element(dof);
                        }

                        // Set the Dirichelet values (hanging node constraints don't need to be distributed every time step, only at output)
                        if (has_Dirichlet_BCs){
                            constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                        }
                        solutionSet[fieldIndex]->update_ghost_values();

                        // Print update to screen
                        if (currentIncrement%userInputs.skip_print_steps==0){
                            sprintf(buffer, "field '%2s' [auxiliary solve]: current solution: %12.6e, current residual:%12.6e\n", \
                            fields[fieldIndex].name.c_str(),				\
                            solutionSet[fieldIndex]->l2_norm(),			\
                            residualSet[fieldIndex]->l2_norm());
                            pcout<<buffer;
                        }

                        // Check to see if this individual variable has converged
                        if (userInputs.var_nonlinear[fieldIndex]){
                            if (userInputs.nonlinear_solver_parameters.getToleranceType(fieldIndex) == ABSOLUTE_SOLUTION_CHANGE){

                                double diff;

                                if (fields[fieldIndex].type == SCALAR){
                                    dU_scalar -= *solutionSet[fieldIndex];
                                    diff = dU_scalar.l2_norm();
                                }
                                else {
                                    dU_vector -= *solutionSet[fieldIndex];
                                    diff = dU_vector.l2_norm();
                                }
                                if (currentIncrement%userInputs.skip_print_steps==0){
                                    pcout << "Relative difference between nonlinear iterations: " << diff << " " << nonlinear_it_index << " " << currentIncrement << std::endl;
                                }

                                if (diff > userInputs.nonlinear_solver_parameters.getToleranceValue(fieldIndex) && nonlinear_it_index < userInputs.nonlinear_solver_parameters.getMaxIterations()){
                                    nonlinear_it_converged = false;
                                }

                            }
                            else {
                                std::cerr << "PRISMS-PF Error: Nonlinear solver tolerance types other than ABSOLUTE_CHANGE have yet to be implemented." << std::endl;
                            }
                        }
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
