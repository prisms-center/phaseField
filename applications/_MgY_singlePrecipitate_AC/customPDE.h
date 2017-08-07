#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
	customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {
		c_dependent_misfit = false;
		for (unsigned int i=0; i<dim; i++){
			for (unsigned int j=0; j<dim; j++){
				if (std::abs(sfts_linear1[i][j])>1.0e-12){
					c_dependent_misfit = true;
				}
			}
		}
	};

private:
	#include "../../include/typeDefs.h"

	const userInputParameters<dim> userInputs;

	// Pure virtual method in MatrixFreePDE
	void residualRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Pure virtual method in MatrixFreePDE
	void residualLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	void solveIncrement();

	// Virtual method in MatrixFreePDE that we override if we need postprocessing
	#ifdef POSTPROCESS_FILE_EXISTS
	void postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
					const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
	#endif

	// Virtual method in MatrixFreePDE that we override if we need nucleation
	#ifdef NUCLEATION_FILE_EXISTS
	double getNucleationProbability(variableValueContainer variable_value, double dV) const;
	#endif

	// ================================================================
	// Methods specific to this subclass
	// ================================================================


	// ================================================================
	// Model constants specific to this subclass
	// ================================================================

	double McV = userInputs.get_model_constant_double("McV");
	double Mn1V = userInputs.get_model_constant_double("Mn1V");
	dealii::Tensor<2,dim> Kn1 = userInputs.get_model_constant_rank_2_tensor("Kn1");
	double W = userInputs.get_model_constant_double("W");
	bool n_dependent_stiffness = userInputs.get_model_constant_bool("n_dependent_stiffness");
	dealii::Tensor<2,dim> sfts_linear1 = userInputs.get_model_constant_rank_2_tensor("sfts_linear1");
	dealii::Tensor<2,dim> sfts_const1 = userInputs.get_model_constant_rank_2_tensor("sfts_const1");
	double A2 = userInputs.get_model_constant_double("A2");
	double A1 = userInputs.get_model_constant_double("A1");
	double A0 = userInputs.get_model_constant_double("A0");
	double B2 = userInputs.get_model_constant_double("B2");
	double B1 = userInputs.get_model_constant_double("B1");
	double B0 = userInputs.get_model_constant_double("B0");

	const static unsigned int CIJ_tensor_size =2*dim-1+dim/3;

	bool c_dependent_misfit;

	double integrated_c_before;

	double dt_modifier;

	// ================================================================

};

//solveIncrement() method for MatrixFreePDE class

//solve each time increment
template <int dim, int degree>
void customPDE<dim,degree>::solveIncrement(){

    //log time
    this->computing_timer.enter_section("matrixFreePDE: solveIncrements");
    Timer time;
    char buffer[200];

    //compute residual vectors
    this->computeRHS();

    //solve for each field
    for(unsigned int fieldIndex=0; fieldIndex<this->fields.size(); fieldIndex++){
        this->currentFieldIndex = fieldIndex; // Used in computeLHS()

        //Parabolic (first order derivatives in time) fields
        if (this->fields[fieldIndex].pdetype==PARABOLIC){

            // Explicit-time step each DOF
            // Takes advantage of knowledge that the length of solutionSet and residualSet is an integer multiple of the length of invM for vector variables
			double nominal_dt;
			if (this->currentIncrement == 1){
				dt_modifier = 0.0;
			}
			else if (this->currentIncrement == 2){
				dt_modifier = 1.0;
			}

			if (fieldIndex == 0 && this->currentIncrement == 1){
				this->computeIntegralMF(integrated_c_before, fieldIndex, this->solutionSet);
			}

            unsigned int invM_size = this->invM.local_size();
            for (unsigned int dof=0; dof<this->solutionSet[fieldIndex]->local_size(); ++dof){
                this->solutionSet[fieldIndex]->local_element(dof)=			\
                this->invM.local_element(dof%invM_size)*this->residualSet[fieldIndex]->local_element(dof);
            }

			if (fieldIndex == 0){
				this->constraintsOtherSet[fieldIndex]->distribute(*(this->solutionSet[fieldIndex]));
				this->constraintsDirichletSet[fieldIndex]->distribute(*(this->solutionSet[fieldIndex]));
	            this->solutionSet[fieldIndex]->update_ghost_values();

				double integrated_c_after;
				this->computeIntegralMF(integrated_c_after, fieldIndex, this->solutionSet);

				double domain_volume = 1.0;
				for (unsigned int i=0; i<dim; i++){
					domain_volume *= userInputs.domain_size[i];
				}

				for (unsigned int dof=0; dof<this->solutionSet[fieldIndex]->local_size(); ++dof){
	                this->solutionSet[fieldIndex]->local_element(dof) -= (integrated_c_after - integrated_c_before)/(domain_volume);
	            }
				//this->pcout << "Before, after, shift:" << integrated_c_before << " " << integrated_c_after << " " << (integrated_c_after - integrated_c_before)/(domain_volume) << std::endl;
			}

            // Set the Dirichelet values (hanging node constraints don't need to be distributed every time step, only at output)
            this->constraintsDirichletSet[fieldIndex]->distribute(*(this->solutionSet[fieldIndex]));
            this->solutionSet[fieldIndex]->update_ghost_values();

            // Print update to screen
            if (this->currentIncrement%userInputs.skip_print_steps==0){
                sprintf(buffer, "field '%2s' [explicit solve]: current solution: %12.6e, current residual:%12.6e\n", \
                this->fields[fieldIndex].name.c_str(),				\
                this->solutionSet[fieldIndex]->l2_norm(),			\
                this->residualSet[fieldIndex]->l2_norm());
                this->pcout<<buffer;
            }
        }
        //Elliptic (time-independent) fields
        else if (this->fields[fieldIndex].pdetype==ELLIPTIC){

            //implicit solve
            //apply Dirichlet BC's
            // Loops through all DoF to which ones have Dirichlet BCs applied, replace the ones that do with the Dirichlet value
            // Is this needed? Why are we applying BCs to the residualSet?
            // This clears the residual where we want to apply Dirichlet BCs, otherwise the solver sees a positive residual
            for (std::map<types::global_dof_index, double>::const_iterator it=this->valuesDirichletSet[fieldIndex]->begin(); it!=this->valuesDirichletSet[fieldIndex]->end(); ++it){
                if (this->residualSet[fieldIndex]->in_local_range(it->first)){
                    (*(this->residualSet[fieldIndex]))(it->first) = 0.0; //it->second; //*jacobianDiagonal(it->first);
                }
            }

            //solver controls
            double tol_value;
            if (userInputs.abs_tol == true){
                tol_value = userInputs.solver_tolerance;
            }
            else {
                tol_value = userInputs.solver_tolerance*this->residualSet[fieldIndex]->l2_norm();
            }

            SolverControl solver_control(userInputs.max_solver_iterations, tol_value);

            // Currently the only allowed solver is SolverCG, the SolverType input variable is a dummy
            SolverCG<vectorType> solver(solver_control);

            //solve
            try{
                if (this->fields[fieldIndex].type == SCALAR){
                    this->dU_scalar=0.0;
                    solver.solve(*this, this->dU_scalar, *(this->residualSet[fieldIndex]), IdentityMatrix(this->solutionSet[fieldIndex]->size()));
                }
                else {
                    this->dU_vector=0.0;
                    solver.solve(*this, this->dU_vector, *(this->residualSet[fieldIndex]), IdentityMatrix(this->solutionSet[fieldIndex]->size()));
                }
            }
            catch (...) {
                this->pcout << "\nWarning: implicit solver did not converge as per set tolerances. consider increasing maxSolverIterations or decreasing solverTolerance.\n";
            }
            if (this->fields[fieldIndex].type == SCALAR){
                *(this->solutionSet[fieldIndex])+=this->dU_scalar;
            }
            else {
                *(this->solutionSet[fieldIndex])+=this->dU_vector;
            }

            if (this->currentIncrement%userInputs.skip_print_steps==0){
                double dU_norm;
                if (this->fields[fieldIndex].type == SCALAR){
                    dU_norm = this->dU_scalar.l2_norm();
                }
                else {
                    dU_norm = this->dU_vector.l2_norm();
                }
                sprintf(buffer, "field '%2s' [implicit solve]: initial residual:%12.6e, current residual:%12.6e, nsteps:%u, tolerance criterion:%12.6e, solution: %12.6e, dU: %12.6e\n", \
                this->fields[fieldIndex].name.c_str(),			\
                this->residualSet[fieldIndex]->l2_norm(),			\
                solver_control.last_value(),				\
                solver_control.last_step(), solver_control.tolerance(), this->solutionSet[fieldIndex]->l2_norm(), dU_norm);
                this->pcout<<buffer;
            }

        }

        //Hyperbolic (second order derivatives in time) fields and general
        //non-linear PDE types not yet implemented
        else{
            this->pcout << "matrixFreePDE.h: unknown field pdetype\n";
            exit(-1);
        }

        //check if solution is nan
        if (!numbers::is_finite(this->solutionSet[fieldIndex]->l2_norm())){
            sprintf(buffer, "ERROR: field '%s' solution is NAN. exiting.\n\n",
            this->fields[fieldIndex].name.c_str());
            this->pcout<<buffer;
            exit(-1);
        }
    }
    if (this->currentIncrement%userInputs.skip_print_steps==0){
        this->pcout << "wall time: " << time.wall_time() << "s\n";
    }
    //log time
    this->computing_timer.exit_section("matrixFreePDE: solveIncrements");

}
